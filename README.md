# Efficiently Finding Shortest Path on 3D Weighted Terrain Surfaces

## Overview

This project provides the implementation of the algorithm for efficiently finding the shortest path on 3D weighted terrain surface. We refer the readers to our paper for more details.

We compared 18 algorithms as follows:

- algorithm FixSP (baseline)
- algorithm LogSP (baseline)
- algorithm Roug-Ref(NoPrunDijk, FixSP, NoEdgSeqConv, NoEffWeig) (baseline)
- algorithm Roug-Ref(NoPrunDijk, FixSP, NoEdgSeqConv, .) (baseline)
- algorithm Roug-Ref(NoPrunDijk, FixSP, ., NoEffWeig) (baseline)
- algorithm Roug-Ref(NoPrunDijk, FixSP, ., .) (baseline)
- algorithm Roug-Ref(NoPrunDijk, ., NoEdgSeqConv, NoEffWeig) (baseline)
- algorithm Roug-Ref(NoPrunDijk, ., NoEdgSeqConv, .) (baseline)
- algorithm Roug-Ref(NoPrunDijk, ., ., NoEffWeig) (baseline)
- algorithm Roug-Ref(NoPrunDijk, ., ., .) (baseline)
- algorithm Roug-Ref(., FixSP, NoEdgSeqConv, NoEffWeig) (baseline)
- algorithm Roug-Ref(., FixSP, NoEdgSeqConv, .) (baseline)
- algorithm Roug-Ref(., FixSP, ., NoEffWeig) (baseline)
- algorithm Roug-Ref(., FixSP, ., .) (baseline)
- algorithm Roug-Ref(., ., NoEdgSeqConv, NoEffWeig) (baseline)
- algorithm Roug-Ref(., ., NoEdgSeqConv, .) (baseline)
- algorithm Roug-Ref(., ., ., NoEffWeig) (baseline)
- algorithm Roug-Ref (our algorithm)

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Dataset

The dataset are stored in "input/" folder.

The datasets are as follows:

- "BH.off" (with dataset size of 279684)
- "EP.off" (with dataset size of 296576)
- "SB.off" (with dataset size of 2048)
- "CP.off" (with dataset size of 2048)
- "PA.off" (with dataset size of 1136)
- "BH_small.off" (a sub-region of "BH.off" with dataset size of 2738)
- "EP_small.off" (a sub-region of "EP.off" with dataset size of 2738)
- "EP_1000768.off" (generated using "EP.off" with dataset size of 1000768)
- "EP_2002080.off" (generated using "EP.off" with dataset size of 2002080)
- "EP_3001050.off" (generated using "EP.off" with dataset size of 3001050)
- "EP_4003072.off" (generated using "EP.off" with dataset size of 4003072)
- "EP_5004800.off" (generated using "EP.off" with dataset size of 5004800)
- "EP_small_10082.off" (generated using "EP_small.off" with dataset size of 10082)
- "EP_small_20000.off" (generated using "EP_small.off" with dataset size of 20000)
- "EP_small_30258.off" (generated using "EP_small.off" with dataset size of 30258)
- "EP_small_40328.off" (generated using "EP_small.off" with dataset size of 40328)
- "EP_small_50562.off" (generated using "EP_small.off" with dataset size of 50562)

Since the file size for the dataset "EP_1000768.off", "EP_2002080.off", "EP_3001050.off", "EP_4003072.off", and "EP_5004800.off" are too large (i.e., they exceed the maximum file size for Github), please download these five files at https://www.dropbox.com/sh/eqwssovnb17rgos/AAASXI-xAnrYsRTvTWSfnBJ6a?dl=0, and put them back in "input/" folder.

Data Format:

We used the .off format in the experiment. The content of the .off file is as follows:

```
OFF

vertices_num faces_num edges_num

1st_vertex_x_coord 1st_vertex_y_coord 1st_vertex_z_coord

2nd_vertex_x_coord 2nd_vertex_y_coord 2nd_vertex_z_coord

......

last_vertex_x_coord last_vertex_y_coord last_vertex_z_coord

1st_face_1st_vertex_ID 1st_face_2nd_vertex_ID 1st_face_3td_vertex_ID

2nd_face_1st_vertex_ID 2nd_face_2nd_vertex_ID 2nd_face_3td_vertex_ID

......

last_face_1st_vertex_ID last_face_2nd_vertex_ID last_face_3td_vertex_ID
```

## Compile command

```
cd src
g++ -o main main.cpp -std=c++11
```

## Run command

```
./main [terrain_data] [epsilon] [removing_value] [calculate_exact_path] [calculate_FixSP]
```

The meaning for each parameter is as follows:

- [terrain_data]: terrain data file name
- [epsilon]: the epsilon value (0 < epsilon <= 1)
- [removing_value]: the removing value (by default = 2)
- [calculate_exact_path]: whether to calculate exact path, 1 for calculate exact path, 0 for not calculate exact path
- [calculate_FixSP]: whether to calculate FixSP, 1 for calculate FixSP, 0 for not calculate FixSP

By default (even if you set both [calculate_exact_path] and [calculate_FixSP] to 0), the project will run algorithm including LogSP without calculating the error ratio compared with the exact path (thus, the distance error compared with the exact path and the calculated path will be infinity).

When you set [calculate_exact_path] to 0 and [calculate_FixSP] to 1, the project will run all algorithms without calculating the error ratio compared with the exact path.

When you set both [calculate_exact_path] and [calculate_FixSP] to 1, the project will run all algorithms and also calculate the error ratio compared with the exact path.

We use Roug-Ref(NoPrunDijk, FixSP, NoEdgSeqConv, NoEffWeig) and set epsilon = 0.05 and the removing value to be 2 to simulate the exact path. So when you run the terrain data with dataset size more than 2000, we strongly encourage you to set [calculate_exact_path] to 0. Otherwise, it will take a very long time to simulate the exact path (which is not the main purpose of this project). You should simulate the exact path ONLY when you need to calculate the error ratio of the calculated path compared with the exact path in the experiment. In addition, when you run the terrain data with dataset size more than 20000, we strongly encourage you to set [calculate_FixSP] to 0. Otherwise, it will take a very long time to run algorithm FixSP (because algorithm FixSP performs very bad when terrain data is large).

An example:

```
./main EP_small.off 1 2 0 1
```

In this example, EP_small.off is the terrain data file, epsilon is 1, removing value is 2, exact path will not be calculated, and algorithm FixSP will be included (thus, it will run all 18 algorithms).

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```
[dataset] [datasize] [epsilon] [removing_value] [building_time (ms)] [query_time_algo1 (ms)] [query_time_algo2_snell_law (ms)] [query_time_algo2_refined (ms)] [query_time_total (ms)] [memory_usage_algo1 (MB)] [memory_usage_algo2_snell_law (MB)] [memory_usage_algo2_refined (MB)] [memory_usage_total (MB)] [snell_law_iteration_count] [distance_error] [distance] [edge_sequence_size]
```

These information will also be shown in the terminal. 

