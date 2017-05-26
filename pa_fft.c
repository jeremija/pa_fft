/*
 * pa_fft is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * pa_fft is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with pa_fft.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <complex.h>
#include <tgmath.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdbool.h>

#include <pulse/simple.h>
#include <pulse/error.h>

#include <fftw3.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <GL/gl.h>

enum w_type {
    WINDOW_TRIANGLE,
    WINDOW_HANNING,
    WINDOW_HAMMING,
    WINDOW_BLACKMAN,
    WINDOW_BLACKMAN_HARRIS,
    WINDOW_WELCH,
    WINDOW_FLAT,
};

struct pa_fft {
    pthread_t thread;
    bool stop, log_graph, overlap, no_refresh;
    int cont;

    /* Pulse */
    pa_simple *s;
    const char *dev;
    int error;
    pa_sample_spec ss;
    pa_channel_map map;

    /* Buffer */
    float *pa_buf;
    size_t pa_buf_size;
    unsigned int pa_samples;
    double *buffer;
    size_t buffer_size;
    unsigned int buffer_samples;
    float **frame_avg_mag;
    unsigned int size_avg;
    unsigned int frame_avg;

    /* FFT */
    int fft_flags;
    unsigned int start_low;
    enum w_type win_type;
    fftw_complex *output;
    unsigned int output_size;
    unsigned int fft_memb;
    double fft_fund_freq;
    fftw_plan plan;

    /* SDL */
    int width, height;
    SDL_Window *win;
    SDL_Event event;
};

void graph_init(struct pa_fft *c)
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    c->win = SDL_CreateWindow("pa_fft", 0, 0, c->width, c->height,
                              SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    SDL_GL_SetAttribute(SDL_GL_ACCELERATED_VISUAL, 1);
    SDL_GL_CreateContext(c->win);
    SDL_GL_SetSwapInterval(1);

    printf("OpenGL Version %s\n", glGetString(GL_VERSION));

    glViewport(0, 0, c->width, c->height);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    SDL_GL_SwapWindow(c->win);
}

void deinit_fft(struct pa_fft *pa_fft) {
    if (!pa_fft)
        return;

    if (pa_fft->cont != 2) {
        pa_fft->cont = 0;
        sleep(1);
    }

    fftw_destroy_plan(pa_fft->plan);
    fftw_free(pa_fft->output);

    if (pa_fft->cont != 2)
        pa_simple_free(pa_fft->s);
    free(pa_fft->buffer);
    free(pa_fft->pa_buf);

    pthread_join(pa_fft->thread, NULL);

    free(pa_fft);
}

static inline void avg_buf_init(struct pa_fft *pa_fft)
{
    pa_fft->frame_avg_mag = malloc(pa_fft->pa_samples*sizeof(float *));
    for (int i = 0; i < pa_fft->pa_samples; i++)
        pa_fft->frame_avg_mag[i] = calloc(pa_fft->frame_avg, sizeof(float));
}

static inline void weights_init(float *dest, int samples, enum w_type w)
{
    switch(w) {
        case WINDOW_TRIANGLE:
            for (int i = 0; i < samples; i++)
                dest[i] = 1 - 2*fabsf((i - ((samples - 1)/2.0f))/(samples - 1));
            break;
        case WINDOW_HANNING:
            for (int i = 0; i < samples; i++)
                dest[i] = 0.5f*(1 - cos((2*M_PI*i)/(samples - 1)));
            break;
        case WINDOW_HAMMING:
            for (int i = 0; i < samples; i++)
                dest[i] = 0.54 - 0.46*cos((2*M_PI*i)/(samples - 1));
            break;
        case WINDOW_BLACKMAN:
            for (int i = 0; i < samples; i++) {
                const float c1 = cos((2*M_PI*i)/(samples - 1));
                const float c2 = cos((4*M_PI*i)/(samples - 1));
                dest[i] = 0.42659 - 0.49656*c1 + 0.076849*c2;
            }
            break;
        case WINDOW_BLACKMAN_HARRIS:
            for (int i = 0; i < samples; i++) {
                const float c1 = cos((2*M_PI*i)/(samples - 1));
                const float c2 = cos((4*M_PI*i)/(samples - 1));
                const float c3 = cos((6*M_PI*i)/(samples - 1));
                dest[i] = 0.35875 - 0.48829*c1 + 0.14128*c2 - 0.001168*c3;
            }
            break;
        case WINDOW_FLAT:
            for (int i = 0; i < samples; i++)
                dest[i] = 1.0f;
            break;
        case WINDOW_WELCH:
            for (int i = 0; i < samples; i++)
                dest[i] = 1 - pow((i - ((samples - 1)/2.0f))/((samples - 1)/2.0f), 2.0f);
            break;
        default:
            for (int i = 0; i < samples; i++)
                dest[i] = 0.0f;
            break;
    }
    float sum = 0.0f;
    for (int i = 0; i < samples; i++)
        sum += dest[i];
    for (int i = 0; i < samples; i++)
        dest[i] /= sum;
}

static inline void apply_win(double *dest, float *src, float *weights,
                             int samples)
{
    for (int i = 0; i < samples; i++)
        dest[i] = src[i]*weights[i];
}

static inline float frame_average(float mag, float *buf, int avgs, int no_mod)
{
    if (!avgs)
        return mag;
    float val = mag;
    for (int i = 0; i < avgs; i++)
        val += buf[i];
    val /= avgs + 1;
    if (no_mod)
        return val;
    for (int i = avgs - 1; i > 0; i--)
        buf[i] = buf[i-1];
    buf[0] = mag;
    return val;
}

void *pa_fft_thread(void *arg) {
    struct pa_fft *t = (struct pa_fft *)arg;
    float weights[t->buffer_samples];

    graph_init(t);
    avg_buf_init(t);
    weights_init(weights, t->fft_memb, t->win_type);

    while (t->cont) {
        while(SDL_PollEvent(&t->event)) {
            switch(t->event.type) {
                case SDL_QUIT:
                    t->cont = 0;
                    break;
                case SDL_WINDOWEVENT:
                    SDL_GL_GetDrawableSize(t->win, &t->width, &t->height);
                    glViewport(0, 0, t->width, t->height);
                    break;
                default:
                    break;
            }
        }

        /* if (t->overlap) */
        /*     memcpy(&t->pa_buf[0], &t->pa_buf[t->pa_samples], */
        /*            t->pa_samples*sizeof(float)); */

        pa_usec_t lag = pa_simple_get_latency(t->s, &t->error);

        if (pa_simple_read(t->s, &t->pa_buf[0],
            t->pa_buf_size, &t->error) < 0) {
            fprintf(stderr, __FILE__": pa_simple_read() failed: %s\n",
                    pa_strerror(t->error));
            t->cont = 0;
            continue;
        }

        /* apply_win(t->buffer, t->pa_buf, weights, t->buffer_samples); */
        /* fftw_execute(t->plan); */

        /* double freq_low, freq_disp, freq_range, freq_off, mag_max = 0.0f; */
        /* if (t->log_graph) { */
        /*     freq_low = log10((t->start_low*t->fft_fund_freq)/((float)t->ss.rate/2)); */
        /*     freq_disp = 1.0 - log10((t->fft_memb*t->fft_fund_freq)/((float)t->ss.rate/2)); */
        /*     freq_range = (1.0 - freq_disp) - freq_low; */
        /*     freq_off = 0.0f; */
        /* } else { */
        /*     freq_low = (t->start_low*t->fft_fund_freq)/((float)t->ss.rate/2); */
        /*     freq_disp = 1.0 - (t->fft_memb*t->fft_fund_freq)/((float)t->ss.rate/2); */
        /*     freq_range = (1.0 - freq_disp) - freq_low; */
        /*     freq_off = 1.0f; */
        /* } */
        /* for (int i = t->start_low; i < t->fft_memb; i++) { */
        /*     fftw_complex num = t->output[i]; */
        /*     double mag = creal(num)*creal(num) + cimag(num)*cimag(num); */
        /*     mag = log10(mag)/10; */
        /*     mag = frame_average(mag, t->frame_avg_mag[i], t->frame_avg, 1); */
        /*     mag_max = mag > mag_max ? mag : mag_max; */
        /* } */

        if (!t->no_refresh)
            glClear(GL_COLOR_BUFFER_BIT);
        glBegin(GL_LINE_STRIP);
        if ((float)lag/1000000 < 1.0f)
            glColor3f(255.0,255.0,255.0);
        else
            glColor3f(255.0,0.0,0.0);
        /* for (int i = t->start_low; i < t->fft_memb; i++) { */
        /*     double freq; */
        /*     fftw_complex num = t->output[i]; */
        /*     if (t->log_graph) */
        /*         freq = log10((i*t->fft_fund_freq)/((float)t->ss.rate/2)); */
        /*     else */
        /*         freq = (i*t->fft_fund_freq)/((float)t->ss.rate/2); */
        /*     double mag = creal(num)*creal(num) + cimag(num)*cimag(num); */
        /*     mag = log10(mag)/10; */
        /*     mag = frame_average(mag, t->frame_avg_mag[i], t->frame_avg, 0); */
        /*     glVertex2f((freq/freq_range + freq_disp/2)*2 - freq_off, mag + mag_max + 0.5f); */
        /* } */
        for (int i = 0; i < t->pa_samples; i++) {
            float mag = t->pa_buf[i];
            mag = frame_average(mag, t->frame_avg_mag[i], t->frame_avg, 0);
            glVertex2f((float) i / t->pa_samples * 2 - 1, (float) mag * 10);
        }
        glEnd();

        SDL_GL_SwapWindow(t->win);
    }

    SDL_DestroyWindow(t->win);
    SDL_Quit();

    deinit_fft(t);

    return NULL;
}

static inline void init_pulse(struct pa_fft *pa_fft)
{
    /* PA spec */
    fprintf(stderr, "device = %s\n", pa_fft->dev);
    if (!pa_fft->dev) {
        fprintf(stderr, "Warning: no device specified! It's highly possible "
                        "Pulseaudio will attempt to use the microphone!\n");
    }

    pa_fft->ss.format = PA_SAMPLE_FLOAT32LE;
    pa_fft->ss.rate = 44100;
    pa_fft->ss.channels = 1;
    pa_channel_map_init_mono(&pa_fft->map);

    if (!(pa_fft->s = pa_simple_new(NULL, "pa_fft", PA_STREAM_RECORD, pa_fft->dev,
                                    "record", &pa_fft->ss, &pa_fft->map, NULL,
                                    &pa_fft->error))) {
        fprintf(stderr, __FILE__": pa_simple_new() failed: %s\n",
                pa_strerror(pa_fft->error));
        pa_fft->cont = 0;
        return;
    }
}

static inline void init_buffers(struct pa_fft *pa_fft)
{
    if (!pa_fft)
        return;

    /* Pulse buffer */
    pa_fft->pa_samples = pa_fft->buffer_samples/(pa_fft->overlap ? 2 : 1);
    pa_fft->pa_buf_size = sizeof(float)*pa_fft->buffer_samples;
    pa_fft->pa_buf = malloc(pa_fft->pa_buf_size);

    /* Input buffer */
    pa_fft->buffer_size = sizeof(double)*pa_fft->buffer_samples;
    pa_fft->buffer = malloc(pa_fft->buffer_size);

    /* FFTW buffer */
    pa_fft->output_size = sizeof(fftw_complex)*pa_fft->buffer_samples;
    pa_fft->output = fftw_malloc(pa_fft->output_size);
    pa_fft->fft_memb = (pa_fft->buffer_samples/2)+1;
    pa_fft->fft_fund_freq = (double)pa_fft->ss.rate/pa_fft->buffer_samples;
}

static inline void init_fft(struct pa_fft *pa_fft) {
    if (!pa_fft)
        return;

    pa_fft->plan = fftw_plan_dft_r2c_1d(pa_fft->buffer_samples, pa_fft->buffer,
                                        pa_fft->output, pa_fft->fft_flags);
}

static inline void print_help()
{
    fprintf(stderr, "Usage: pa_fft <options>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -l        Linear (instead of logarithmic) graph\n");
    fprintf(stderr, "    -o        Use overlap\n");
    fprintf(stderr, "    -n        Do not clear window on every frame\n");
    fprintf(stderr, "    -h        Print this\n");
    fprintf(stderr, "    -a (int)  Turn averaging over (int) frames on (default = 1 frame)\n");
    fprintf(stderr, "    -c (int)  Cut off low frequency coeffs (default = 1st coef onward)\n");
    fprintf(stderr, "    -s (int)  Amount of samples to use for transform (default = 1024)\n");
    fprintf(stderr, "    -w (str)  Specify windowing function (default = \"hanning\", possible values:\n");
    fprintf(stderr, "                  \"triangle\" \"hanning\" \"hamming\" \"blackman\""
                                      "\"blackman-harris\" \"welch\" \"flat\"\n");
    fprintf(stderr, "    -d (str)  Pulseaudio device\n");
    fprintf(stderr, "                  Specify using the name from \"pacmd list-sources | grep \"name:\"\"\n");
}

int main(int argc, char *argv[])
{
    int c;
    struct pa_fft *ctx = calloc(1, sizeof(struct pa_fft));
    ctx->cont = 1;

    if (argc < 2 || !strcmp(argv[1], "-h")) {
        print_help();
        return 1;
    }

    ctx->width = 1200;
    ctx->height = 500;
    ctx->buffer_samples = 1024;
    ctx->dev = NULL;
    ctx->log_graph = 1;
    ctx->no_refresh = 0;
    ctx->overlap = 0;
    ctx->frame_avg = 2;
    ctx->start_low = 1;
    ctx->win_type = WINDOW_BLACKMAN_HARRIS;
    ctx->fft_flags = FFTW_PATIENT | FFTW_DESTROY_INPUT;

    const char *opt_str = "lonhd:a:c:s:w:";
    while ((c = getopt (argc, argv, opt_str)) != -1) {
        switch (c) {
            case 'l':
                ctx->log_graph = 0;
                break;
            case 'o':
                ctx->overlap = 1;
                break;
            case 'n':
                ctx->no_refresh = 1;
                break;
            case 'd':
                ctx->dev = strdup(optarg);
                break;
            case 'h':
                print_help();
                return 0;
                break;
            case 'a':
                sscanf(optarg, "%u", &ctx->frame_avg);
                break;
            case 'c':
                sscanf(optarg, "%u", &ctx->start_low);
                break;
            case 's':
                sscanf(optarg, "%u", &ctx->buffer_samples);
                break;
            case 'w':
                if (!strcmp(optarg, "triangle"))
                    ctx->win_type = WINDOW_TRIANGLE;
                else if (!strcmp(optarg, "hanning"))
                    ctx->win_type = WINDOW_HANNING;
                else if (!strcmp(optarg, "hamming"))
                    ctx->win_type = WINDOW_HAMMING;
                else if (!strcmp(optarg, "blackman"))
                    ctx->win_type = WINDOW_BLACKMAN;
                else if (!strcmp(optarg, "blackman-harris"))
                    ctx->win_type = WINDOW_BLACKMAN_HARRIS;
                else if (!strcmp(optarg, "welch"))
                    ctx->win_type = WINDOW_WELCH;
                else if (!strcmp(optarg, "flat"))
                    ctx->win_type = WINDOW_FLAT;
                else
                    fprintf(stderr, "Unknown window \"%s\"\n", optarg);
                break;
            case '?':
                for (int i = 0; i < sizeof(opt_str); i++) {
                    if (opt_str[i] == optopt) {
                        fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                    }
                }
            default:
                abort();
        }
    }

    init_pulse(ctx);
    init_buffers(ctx);
    init_fft(ctx);

    if (!ctx->cont) {
        ctx->cont = 2;
        deinit_fft(ctx);
        return 1;
    }

    pthread_create(&ctx->thread, NULL, pa_fft_thread, ctx);

    while(ctx->cont) {
        sleep(1);
    }

    return 0;
}
