# run15_pid_in_jets

## Read pico DSTs and save tracks/jets into trees
1. Edit `src/tracks.cxx` or `src/jets.cxx`
2. Enter container and compile with `source runimage.sh`, `source setup.sh` and `make`
3. Submit batch jobs with `star-submit submit.xml`

## Analyze output trees
`root -l analysis_code/analyze_track_tree.C` or `root -l analysis_code/analyze_jet_tree.C`
