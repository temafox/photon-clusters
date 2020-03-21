# photon-clusters
Code for my third-year coursework

Structure of the project:
* This README file for documentation
* Data generation
	* src/data-gen/makeData.cpp - generate a tree of model signal
	* src/data-gen/Makefile - via **make**, build the engine and make it create trees
	* data/ - a directory for data generation and access to data
* Algorithm training and testing
	* src/algorithm/example/ - an example of algorithm training and testing
		* tmvaClassification.cpp - training and testing of the TMVA classifiers, placing the output into *data/tmva.root*
		* tmvaClassificationApplication.cpp - an analysis module to use the trained classifiers
* Startup scripts
	* scripts/generate - generate a training tree and a testing tree
	* scripts/classify - given a training tree and a testing tree in *data/*, launch the *tmvaClassification.cpp* macro
	* scripts/app - given *data/tmva.root*, launch the *tmvaClassificationApplication.cpp* module
