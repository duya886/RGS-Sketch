# RGS-Sketch

RGS-Sketch is a novel invertible sketch that simultaneously achieves high memory efficiency, supports online detection, and offers mergeability. It is submitted to VLDB 2025. This is the repository for the source codes.

## Abstract
Super spreader detection in high-speed data streams is crucial for numerous applications. Among various methods, invertible methods are getting more attention since they can recover the IDs of super spreaders just from their own data structures. Nevertheless, there remains a lack of an invertible approach that can concurrently adapt to skewed data streams, support online detection, and enable merging data from different measurement points/periods. This paper proposes RGS-Sketch, a novel invertible sketch designed to address these challenges. The core of RGS-Sketch lies in a new mergeable memory sharing design called register group sharing. This design organizes registers (small-sized counters) into groups as basic memory-sharing units, accommodating the high skewness of real-world data streams and offering high memory efficiency. Besides, it enables online cardinality estimation through the real-time acquisition of a group's state. To enhance detection accuracy further, we propose a limited register update strategy. It blocks flows with small cardinalities from updating registers, thereby reducing memory overhead and estimation noises. Extensive experimental results based on four real-world datasets show that RGS-Sketch significantly outperforms the most accurate baselines in accuracy while maintaining a high throughput. Specifically, it improves the F1 Scores by up to 0.643 for measurements at a single point/period and up to 0.472 for measurements across multiple points/periods.

## About this repo
The codes of RGS-Sketch and the baselines are compiled using gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04). For RGS-Sketch, FreeRS, SpreadSketch, ExtendedSketch+, FlowFight, the compiling command is as follows (“***“ represents the name of the main algorithm code file):
```shell
g++ ***.cpp MurmurHash3.cpp -m64 -O3
```

For ExtendedSketch, the compiling command is as follows:
```shell
g++ ExtendedSketch.cpp MurmurHash3.cpp mini-gmp.c -m64 -O3
```



On win10 with MinGW, the compiling commands may need to add " -static-libstdc++ -static-libgcc".

