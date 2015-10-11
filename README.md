## pa_fft

--------------

##A pulseaudio spectrum analyzer

Requires the OpenGL, SDL2, pthreads and FFTW libraries.

To compile just type `make` and be sure to have all libraries required.

To get information on how to use it type `./pa_fft` or add the `-h` option.

To list all pulseaudio devices use `pacmd list-sources | grep 'name:'`.
Remember to use the name given inside the <> brackets.

#Features:
- logarithmic display
- averaging
- device selection
- variable sample size
- low frequency cutoff
- resizable window

##Examples:
- ./pa_fft -d alsa_output.usb-C-Media_Electronics_Inc._USB_Advanced_Audio_Device-00.analog-stereo.monitor
- ./pa_fft -d alsa_output.pci-0000_00_1b.0.analog-stereo.monitor -l -a 10

You probably want to use devices ending with `.monitor` to analyze what's currently playing.
The FFT line output will turn red in case the latency gets over 1 second, but it's usually below 0.1 seconds.
