# Warnings to turn off
warning('off', 'Octave:shadowed-function');
warning('off', 'Octave:legacy-function');

# Packages to install
pkg install -forge statistics statistics-resampling

# Packages to load
pkg load statistics statistics-resampling

# Warnings to turn back on
warning('on', 'Octave:shadowed-function');

