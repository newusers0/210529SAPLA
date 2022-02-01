<!-- An Indexable Time Series Dimensionality Reduction Method for Max Deviation Reduction and Similarity Search -->
---This repository contains VS2019 C++ implementation of the [SAPLA](https://sites.google.com/view/sapla0/home).

## Introduction
```
This is implementation of SAPLA, APLA, APCA, PLA, PAA, CHEBY, PAALM, SAX, R-tree, updated R-tree, KNN.
```
* SAPLA algorithm is implemented in [CAPLA](https://github.com/newusers0/210529SAPLA/blob/master/CAPLA.h).
* APLA algorithm is implemented in [APLA_ICDE07](https://github.com/newusers0/210529SAPLA/blob/master/APLA_ICDE07.h).
* APCA algorithm is implemented in [CAPCA](https://github.com/newusers0/210529SAPLA/blob/master/CAPCA.h).
* PLA algorithm is implemented in [CPLA](https://github.com/newusers0/210529SAPLA/blob/master/CPLA.h).
* PAA algorithm is implemented in [CAPCA](https://github.com/newusers0/210529SAPLA/blob/master/CAPCA.h).
* CHEBY algorithm is implemented in [CCHEBYSHEV](https://github.com/newusers0/210529SAPLA/blob/master/CCHEBYSHEV.h).
* PAALM algorithm is implemented in [CAPCA](https://github.com/newusers0/210529SAPLA/blob/master/CAPCA.h).
* SAX algorithm is implemented in [saxquantizer](https://github.com/newusers0/210529SAPLA/blob/master/lib/saxquantizer.hpp) from original [code](https://github.com/melsabagh/sax).
* R-tree algorithm is implemented in [RTree](https://github.com/newusers0/210529SAPLA/blob/master/lib/RTree.h) from original [code](https://superliminal.com/sources/).
* DBCH-tree algorithm is implemented in [DBCHTree](https://github.com/newusers0/210529SAPLA/blob/master/RTree_partition.h).
* KNN algorithm is implemented in [MULTI](https://github.com/newusers0/210529SAPLA/blob/master/MULTI_DIMENSION.h).
---

## Data

UCR archive 2018: Raw data can be downloaded from [UCRArchive2018](https://www.cs.ucr.edu/~eamonn/time_series_data_2018/).

---
## Dependencies

* C++ compiler with C++11 support
* Boost Math

---
## Configuration
1. Download [Boost](https://www.boost.org/) and include in VS2019
2. Download [VLD](https://github.com/KindDragon/vld) and install it in VS2019.
3. If VS2019 has C4996 warning, turn off [it](https://docs.microsoft.com/en-us/cpp/error-messages/compiler-warnings/compiler-warning-level-3-c4996?f1url=https%3A%2F%2Fmsdn.microsoft.com%2Fquery%2Fdev15.query%3FappId%3DDev15IDEF1%26l%3DEN-US%26k%3Dk(C4996)%26rd%3Dtrue&view=msvc-160&viewFallbackFrom=vs-2019).

---
## Running
```
cd C:
mkdir Code
```
- Move **UCRArchive_2018** to Code folder
- Switch to **Release** mode in VS2019.
- Press **F5**.

---
## Licence
```
MIT
```

---

