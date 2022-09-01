# Finding Shortest Path on 3D Weighted Terrain Surface using Divide-and-Conquer and Effective Weight

## Overview

This project provides the implementation of the algorithm for finding the shortest path on 3D weighted terrain surface using divide-and-conquer and effective weight.

Our algorithm DLSP-EWSL, and the baseline algorithms, i.e., algorithm FSP, algorithm LSP and algorithm FSP-BSSL, are studied in the experiments. But the calculated weighted shortest path of algorithm FSP and algorithm LSP do not follow Snell's law, so they are not our main focus. In order to conduct the ablation study, i.e., show the superior performance of our algorithms in both edge sequence finding step and edge sequence based weighted shortest path finding step, we also interchanged two steps of algorithm DLSP-EWSL and algorithm FSP-BSSL. That is, we also studied algorithm FSP-EWSL and algorithm DLSP-BSSL in the experiments. In total, we compared six algorithms, i.e., algorithm FSP, algorithm LSP, algorithm FSP-BSSL, algorithm FSP-EWSL, algorithm DLSP-BSSL and algorithm DLSP-EWSL. We refer the readers to our paper for more details.

In total, we compared six algorithms as follows:

- algorithm FSP (baseline)
- algorithm LSP (baseline)
- algorithm FSP-BSSL (baseline)
- algorithm FSP-EWSL (variation)
- algorithm DLSP-BSSL (variation)
- algorithm DLSP-EWSL (our algorithm)

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Dataset

The dataset are stored in "input/" folder.

The datasets are as follows:

- "BH.off" (with dataset size of 279684)
- "EP.off" (with dataset size of 296576)
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
---
```
cd src
g++ -o main main.cpp -std=c++11
```

## Run command

```
./main [terrain_data] [epsilon_SP] [epsilon_SL] [calculate_exact_path] [calculate_FSP]
```

The meaning for each parameter is as follows:

- [terrain_data]: terrain data file name
- [epsilon_SP]: the epsilon value for both DLSP and LSP (0 < epsilon_SP <= 1)
- [epsilon_SL]: the epsilon value for both EWSL and BSSL (0 < epsilon_SL <= 1)
- [calculate_exact_path]: whether to calculate exact path, 1 for calculate exact path, 0 for not calculate exact path
- [calculate_FSP]: whether to calculate FSP, 1 for calculate FSP, 0 for not calculate FSP

By default (even if you set both [calculate_exact_path] and [calculate_FSP] to 0), the project will run algorithm LSP, algorithm DLSP-BSSL and algorithm DLSP-EWSL without calculating the error ratio compared with the exact path (thus, the distance error compared with the exact path and the calculated path will be infinity).

When you set [calculate_exact_path] to 0 and [calculate_FSP] to 1, the project will run algorithm FSP, algorithm LSP, algorithm FSP-BSSL, algorithm FSP-EWSL, algorithm DLSP-BSSL and algorithm DLSP-EWSL without calculating the error ratio compared with the exact path (thus, the distance error compared with the exact path and the calculated path will be infinity).

When you set both [calculate_exact_path] and [calculate_FSP] to 1, the project will run algorithm FSP, algorithm LSP, algorithm FSP-BSSL, algorithm FSP-EWSL, algorithm DLSP-BSSL and algorithm DLSP-EWSL and also calculate the error ratio compared with the exact path.

We use FSP-BSSL and set epsilon = 0.05 (where epsilon = max(epsilon_SP, epsilon_SL) and we set epsilon_SP = epsilon_SL = 0.05) to simulate the exact path. So when you run the terrain data with dataset size more than 2000, we strongly encourage you to set [calculate_exact_path] to 0. Otherwise, it will take a very long time to simulate the exact path (which is not the main purpose of this project). You should simulate the exact path ONLY when you need to calculate the error ratio of the calculated path compared with the exact path in the experiment. In addition, when you run the terrain data with dataset size more than 20000, we strongly encourage you to set [calculate_FSP] to 0. Otherwise, it will take a very long time to run algorithm  FSP (because algorithm FSP performs very bad when terrain data is large).

An example:

```
./main EP_small.off 0.5 0.5 0 1
```

In this example, EP_small.off is the terrain data file, epsilon_SP is 0.5, epsilon_SL is 0.5, exact path will not be calculated (and the distance error compared with the exact path and the calculated path will be infinity), and algorithm FSP will be included (thus, it will run six algorithms, i.e., algorithm FSP, algorithm LSP, algorithm FSP-BSSL, algorithm FSP-EWSL, algorithm DLSP-BSSL and algorithm DLSP-EWSL).

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```
[dataset] [datasize] [epsilon] [epsilon_SP] [epsilon_SL] [building_time (ms)] [query_time_step1 (ms)] [query_time_step2 (ms)] [query_time_total (ms)] [memory_usage_step1 (MB)] [memory_usage_step2 (MB)] [memory_usage_total (MB)] [snell_law_iteration_count] [distance_error] [distance] [edge_sequence_size]
```

These information will also be shown in the terminal. 

