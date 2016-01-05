# build
```sh:qiita.sh
cd build
cmake .. && gmake
```

# run
```sh:qiita.sh
bin/run data/param/parameter_file
```

# parameter file
shell script
* size: system size (size x size)
* N:    step
* eta:  p = C * \nabla \phi ^ \eta
