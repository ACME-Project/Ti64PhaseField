# Ti64PhaseField
Phase field codes for modeling phase transformations in Ti-6Al-4V

# Batch 1 

**Goal** : To implement a vanilla toth model in 3D, without MMSP. Grid size in the source code is 40x40x40

**Files** : NonMMSP_3D.cpp , Tensor.hpp, Eigenstrain.hpp, GradECoeff.hpp

**Compilation** : g++ NonMMSP_3D.cpp -o output

**Purpose** : Models the growth of a spherical nucleus in one of 12 unique directions. The particular variant can be selected in 
          NonMMSP_3D.cpp by choosing the apprpriate epsi'x' variable.
