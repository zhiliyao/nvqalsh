# Introduction
Source code for paper:

Yao Z, Zhang J, Feng J. NV-QALSH: An NVM-Optimized Implementation of Query-Aware Locality-Sensitive Hashing[C]//International Conference on Database and Expert Systems Applications. Springer, Cham, 2021: 58-69.

# How to run
1. Mount an NVM file system on `/home/pmem/test`, or change the path of NVM file system in `main.cc`(change content at `main.cc:243`, `main.cc:247`, `main.cc:251`).

2. Run `exp.sh`, you can find more details in that shell script.

# Notice
NV-QALSH is developed on the basis of [QALSH_Mem](https://github.com/HuangQiang/QALSH_Mem) and [LB-TREE](https://github.com/schencoding/lbtree). They are source codes of the following papers respectively:

Qiang Huang, Jianlin Feng, Yikai Zhang, Qiong Fang, and Wilfred Ng. 2015. Query-aware locality-sensitive hashing for approximate nearest neighbor search. Proc. VLDB Endow. 9, 1 (September 2015), 1–12.

Jihang Liu, Shimin Chen. LB+-Trees: Optimizing Persistent Index Performance on 3DXPoint Memory. PVLDB 13(7): 1078-1090.

We mainly modify the lbtree.h and lbtree.cc to optimize the NVM-based B+-Trees for QALSH, and the source codes of NV-QALSH are mainly the qalsh_tree.h and qalsh_tree.cc