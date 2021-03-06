# photon-clusters
Code for my third-year coursework and diploma

Structure of the project:
* This README file for documentation
* Cluster division
	* src/div/selectPhotons.C - clone a tree and select only those events where a photon appears in the decay of a pion
* Data generation
	* src/data-gen/makeData.cpp - generate a tree of model signal
	* data/ - a directory for data generation and access to data
* Algorithm training and testing
	* src/example/ - an example of algorithm training and testing
		* tmvaClassification.cpp - training and testing of the TMVA classifiers, placing the output into *data/tmva.root*
		* tmvaClassificationApplication.cpp - an analysis module to use the trained classifiers
* Startup scripts
	* scripts/generate *number_of_events* - generate a training tree and a testing tree
	* scripts/train\_classification *[train_test_file]* - given a training tree and a testing tree in *train_test_file* (by default *data/data.root*), launch the *tmvaClassification.cpp* macro, which writes into *data/tmva.root*; also, it writes the best cut value in *data/threshold.Float_t*
	* scripts/show\_classification *[tmva_file]* - show classificators in *tmva_file* if given, otherwise in *data/tmva.root*
	* scripts/apply\_classification *[data_file]* - given *data_file* (by default *data/data.root*) and *data/threshold.Float_t*, launch the *tmvaClassificationApplication.cpp* module, which writes into *data/tmva_app.root*
