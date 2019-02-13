coarse grid
Parametrisierung PARXC PARYC TMAXC
     4  9   0   4     1       NEL NVT NMT NVE NBCT
DCORVG
0 0
0.5 0
1 0
2 0
2.33 0
2.5 0
2.67 0
3 0
1 2
KVERT
 1 2 9 8
 2 3 4 9
 9 4 5 6
6 7 8 9
KNPR
1  1  1  1  1  1  1  1  0
KMM
1 8


#NEL - number of elements in initial mesh
#NVT - number of vertices in the initial domain
#NMT - 
#NVE - number of vertices per element
#NBCT - number of boundary parts
#KVERT - vertices numbers
#KNPR - boundary or interior (1,0)
#KMM - edge which is not mentioned in KVERT