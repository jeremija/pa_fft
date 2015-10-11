## pa_fft

--------------

##A pulseaudio spectrum analyzer

Requires the OpenGL, SDL2 and FFTW libraries.

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
