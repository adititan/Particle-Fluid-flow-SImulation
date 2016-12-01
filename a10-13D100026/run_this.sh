#!/bin/bash

#change hdx = 0.5
python final_assign.py --kernel=CubicSpline --nx=25 --openmp -d ./cubic_1
python final_assign.py --kernel=CubicSpline --nx=50 --openmp -d ./cubic_2
python final_assign.py --kernel=CubicSpline --nx=100 --openmp -d ./cubic_3
python final_assign.py --kernel=QuinticSpline --nx=25 --openmp -d ./quad_1
python final_assign.py --kernel=QuinticSpline --nx=50 --openmp -d ./quad_2
python final_assign.py --kernel=QuinticSpline --nx=100 --openmp -d ./quad_3
python final_assign.py --kernel=Gaussian --nx=25 --openmp -d ./gauss_1
python final_assign.py --kernel=Gaussian --nx=50 --openmp -d ./gauss_2
python final_assign.py --kernel=Gaussian --nx=100 --openmp -d ./gauss_3
python final_assign.py --kernel=WendlandQuintic --nx=25 --openmp -d ./wenl_1
python final_assign.py --kernel=WendlandQuintic --nx=50 --openmp -d ./wenl_2
python final_assign.py --kernel=WendlandQuintic --nx=100 --openmp -d ./wenl_3

#change hdx = 1.0
'python final_assign.py --kernel=CubicSpline --nx=25 --openmp -d ./cubic_4
python final_assign.py --kernel=CubicSpline --nx=50 --openmp -d ./cubic_5
python final_assign.py --kernel=CubicSpline --nx=100 --openmp -d ./cubic_6
python final_assign.py --kernel=QuinticSpline --nx=25 --openmp -d ./quad_4
python final_assign.py --kernel=QuinticSpline --nx=50 --openmp -d ./quad_5
python final_assign.py --kernel=QuinticSpline --nx=100 --openmp -d ./quad_6
python final_assign.py --kernel=Gaussian --nx=25 --openmp -d ./gauss_4
python final_assign.py --kernel=Gaussian --nx=50 --openmp -d ./gauss_5
python final_assign.py --kernel=Gaussian --nx=100 --openmp -d ./gauss_6
python final_assign.py --kernel=WendlandQuintic --nx=25 --openmp -d ./wenl_4
python final_assign.py --kernel=WendlandQuintic --nx=50 --openmp -d ./wenl_5
python final_assign.py --kernel=WendlandQuintic --nx=100 --openmp -d ./wenl_6

#change hdx = 2.0
'
'
python final_assign.py --kernel=CubicSpline --nx=25 --openmp -d ./cubic_7
python final_assign.py --kernel=CubicSpline --nx=50 --openmp -d ./cubic_8
python final_assign.py --kernel=CubicSpline --nx=100 --openmp -d ./cubic_9
python final_assign.py --kernel=QuinticSpline --nx=25 --openmp -d ./quad_7
python final_assign.py --kernel=QuinticSpline --nx=50 --openmp -d ./quad_8
python final_assign.py --kernel=QuinticSpline --nx=100 --openmp -d ./quad_9
python final_assign.py --kernel=Gaussian --nx=25 --openmp -d ./gauss_7
python final_assign.py --kernel=Gaussian --nx=50 --openmp -d ./gauss_8
python final_assign.py --kernel=Gaussian --nx=100 --openmp -d ./gauss_9
python final_assign.py --kernel=WendlandQuintic --nx=25 --openmp -d ./wenl_7
python final_assign.py --kernel=WendlandQuintic --nx=50 --openmp -d ./wenl_8
python final_assign.py --kernel=WendlandQuintic --nx=100 --openmp -d ./wenl_9


'
