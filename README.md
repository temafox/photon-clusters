# photon-clusters
Code for my third-year coursework

Structure of the project:
* This README file for documentation
* Data generation
	* src/data-gen/makeData.cpp - generate a tree of model signal
	* data/ - a directory for data generation and access to data
* Algorithm training and testing
	* src/algorithm/example/ - an example of algorithm training and testing
		* tmvaClassification.cpp - training and testing of the TMVA classifiers, placing the output into *data/tmva.root*
		* tmvaClassificationApplication.cpp - an analysis module to use the trained classifiers
* Startup scripts
	* generate *number_of_events* - generate a training tree and a testing tree
	* classify - given a training tree and a testing tree in *data/data.root*, launch the *tmvaClassification.cpp* macro, which writes into *data/tmva.root*
	* app - given *data/data.root*, launch the *tmvaClassificationApplication.cpp* module, which writes into *data/tmva_app.root*
