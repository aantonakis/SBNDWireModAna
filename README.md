# SBNDWireModAna

## Make WireMod Histograms on the Grid

- 1. Producing histograms with various dimensions on the grid can be done with ``/grid/submit_multi_dim_tracks.py``

This grid submission script is configured with ``/grid/bin/grid_executable_multi_dim_tracks.sh``

- Open the grid executable script to see all available options

- Calibration Options: SCE, YZ Uniformity, E-Lifetime and Calibration Constant

- Various selection options are available (types of cuts to apply)

- The user provides a vector of integers to provide the dimensions to store (10 dimensions available)

- Outputs two types of histograms: One stores the number of hits and the other stores the number of tracks

Example Useage\

```
$ python submit_multi_dim_tracks.py -l file_list.list -o hist_dir_name -nfile 100 -ngrid 10
```

## Project and Merge Histograms on the Grid

- 2. Producing merged projections from step 1 can be done with ``/grid/submit_merge_hists.py``

Similarly, this grid submission script is configured with ``/grid/bin/grid_executable_merge_hists.sh``

- Takes input files from step 1 and projects them down to lower dimensional histograms

- The number of grid jobs is equal to the number of merged outputs from the input file list


Example Useage\

```
$ python submit_merge_hists.py -l file_list.list -o hist_dir_name -nfile 100 -ngrid 10
```


