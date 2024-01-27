# Efficiently Finding Shortest Path on 3D Weighted Terrain Surfaces for Moving Objects

## Overview

This project provides the implementation of the algorithm for efficiently finding the shortest path on 3D weighted terrain surface for moving objects. We refer the readers to our paper for more details.

We compared 13 algorithms as follows:

- algorithm EdgSeq (baseline)
- algorithm FixSP (baseline)
- algorithm FixSP-NoWei-Adp (baseline)
- algorithm LogSP-Adp (baseline)
- algorithm Roug-Ref-Naive1 (baseline)
- algorithm Roug-Ref-NoEffSP (baseline)
- algorithm Roug-Ref-Naive2 (baseline)
- algorithm Roug-Ref-NoEdgSeqConv (baseline)
- algorithm Roug-Ref-NoEffWeig (baseline)
- algorithm Roug-Ref-NoPrunDijk (baseline)
- algorithm LogSP (baseline)
- algorithm Roug (baseline)
- algorithm Roug-Ref (our algorithm)

For algorithm LogSP and Roug, they have distance errors larger than (1+epsilon). In order to compare them with algorithm Roug-Ref, given epsilon, we first calculate a path with distance d. Then, for algorithm LogSP and Roug, we finetune their input error to make their calculated path also with distance d, and measure their performance. 

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Dataset

The dataset are stored in "input/" folder.

The datasets are as follows:

- "BH.off" (with dataset size of 279684)
- "EP.off" (with dataset size of 296576)
- "SB.off" (with dataset size of 2048)
- "VS.off" (with dataset size of 2048)
- "PA.off" (with dataset size of 1136)
- "BH_small.off" (a sub-region of "BH.off" with dataset size of 2738)
- "EP_small.off" (a sub-region of "BH.off" with dataset size of 2738)
- "BH_1001670.off" (generated using "BH.off" with dataset size of 1001670)
- "BH_2003274.off" (generated using "BH.off" with dataset size of 2003274)
- "BH_3002982.off" (generated using "BH.off" with dataset size of 3002982)
- "BH_4004364.off" (generated using "BH.off" with dataset size of 4004364)
- "BH_5002604.off" (generated using "BH.off" with dataset size of 5002604)
- "BH_small_10082.off" (generated using "BH_small.off" with dataset size of 10082)
- "BH_small_20000.off" (generated using "BH_small.off" with dataset size of 20000)
- "BH_small_30258.off" (generated using "BH_small.off" with dataset size of 30258)
- "BH_small_40328.off" (generated using "BH_small.off" with dataset size of 40328)
- "BH_small_50562.off" (generated using "BH_small.off" with dataset size of 50562)
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

Since the file size for the dataset "BH_1001670.off", "BH_2003274.off", "BH_3002982.off", "BH_4004364.off", "BH_5002604.off", "EP_1000768.off", "EP_2002080.off", "EP_3001050.off", "EP_4003072.off", and "EP_5004800.off" are too large (i.e., they exceed the maximum file size for Github), please download these five files at https://www.dropbox.com/sh/eqwssovnb17rgos/AAASXI-xAnrYsRTvTWSfnBJ6a?dl=0, and put them back in "input/" folder.

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
./main [terrain_data] [epsilon] [removing_value] [calculate_exact_path] [calculate_FixSP_comp_and_LogSP_and_Roug]
```

The meaning for each parameter is as follows:

- [terrain_data]: terrain data file name
- [epsilon]: the epsilon value (0 < epsilon <= 1)
- [removing_value]: the removing value (by default = 2)
- [calculate_exact_path]: whether to calculate exact path, 1 for calculate exact path, 0 for not calculate exact path (in this case, the distance error compared with the exact path and the calculated path will be infinity)
- [calculate_FixSP_comp_and_LogSP_and_Roug]: whether to calculate algorithms involves FixSP and algorithm LogSP (with same calculated distance of Roug-Ref) and algorithm Roug (with same calculated distance of Roug-Ref), 1 for calculate, 0 for not calculate

We use algorithm EdgSeq (which involves FixSP) and set epsilon = 0.05 to simulate the exact path. So when you run the terrain data with dataset size more than 2000, we strongly encourage you to set [calculate_exact_path] to 0. Otherwise, it will take a very long time to simulate the exact path (which is not the main purpose of this project). You should simulate the exact path ONLY when you need to calculate the error ratio of the calculated path compared with the exact path in the experiment. In addition, when you run the terrain data with dataset size more than 20000, we strongly encourage you to set [calculate_FixSP_comp_and_LogSP_and_Roug] to 0. Otherwise, it will take a very long time to run algorithm EdgSeq, FixSP and FixSP-NoWei-Adp, LogSP (with same calculated distance of Roug-Ref) and Roug (with same calculated distance of Roug-Ref), because they perform very bad when terrain data is large. 

An example:

```
./main EP_small.off 1 2 0 1
```

In this example, EP_small.off is the terrain data file, epsilon is 1, removing value is 2, exact path will not be calculated, and algorithm EdgSeq, FixSP, FixSP-NoWei-Adp, Roug-Ref-Naive1, Roug-Ref-NoEffSP, LogSP (with same calculated distance of Roug-Ref) and Roug (with same calculated distance of Roug-Ref) will be included. Thus, it will run all 13 algorithms.

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```
[dataset] [datasize] [epsilon] [removing_value] [building_time (ms)] [query_time_algo1 (ms)] [query_time_algo2_snell_law (ms)] [query_time_algo2_refined (ms)] [query_time_total (ms)] [memory_usage_algo1 (MB)] [memory_usage_algo2_snell_law (MB)] [memory_usage_algo2_refined (MB)] [memory_usage_total (MB)] [chances_of_using_error_guaranteed_path_refinement_step] [average_number_of_Steiner_points_per_edge_algo1] [average_number_of_Steiner_points_per_edge_algo2] [total_average_number_of_Steiner_points_per_edge] [snell_law_iteration_count] [distance_error] [distance] [edge_sequence_size]
```

These information will also be shown in the terminal. 

