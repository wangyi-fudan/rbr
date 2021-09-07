all: rbr_model rbr_predict
rbr_model: rbr_model.cpp
	g++ rbr_model.cpp -o rbr_model -O3 -Wall -fopenmp -lgsl -lgslcblas -static -s
rbr_predict: rbr_predict.cpp
	g++ rbr_predict.cpp -o rbr_predict -O3 -Wall -fopenmp -lgsl -lgslcblas -static -s

