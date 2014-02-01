import numpy as np
import copy, sys, os
import shutil
import math
import yaml
from sklearn.externals import joblib
from sklearn.svm import SVC
from sklearn import tree
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import LogisticRegression

from sklearn.cross_validation import KFold
from collections import Counter
from scipy.stats import mannwhitneyu
import json



class classifiers():

    classifierName = ""
    classifierParams = ""

    def __init__(self, classifierName="", classifierParams="default"):
        print "Initializing classifier class..."

        if len(classifierName) > 0: #set the classifier type in the class variable
            self.setClassifierName(classifierName)
            self.setClassifierParams(classifierParams)
        else:
            print "No classifier set yet"

    def performTrainTest(self, trainFeatures, trainClassLabels, evalFeatures, classifierName = "", classifierParams = "default"):

        if len(classifierName) > 0:
            self.classifierName = classifierName
        else:
            if len(self.classifierName)<0:
                print "Please set a valid classifier while initialization of the class, or pass into this function"
                return "Classifier not set"

        if isinstance(classifierParams,dict):
            self.classifierParams = classifierParams

        if self.classifierName == "svm":
            return self.svm(trainFeatures, trainClassLabels, evalFeatures,self.classifierParams)
        elif self.classifierName == "tree":
            return self.tree(trainFeatures, trainClassLabels, evalFeatures,self.classifierParams)
        elif self.classifierName =="mYkNN":
            return self.mYkNN(trainFeatures, trainClassLabels, evalFeatures,self.classifierParams)
        elif self.classifierName =="kNN":
            return self.kNN(trainFeatures, trainClassLabels, evalFeatures,self.classifierParams)
        elif self.classifierName =="nbMulti":
            return self.nbMulti(trainFeatures, trainClassLabels, evalFeatures,self.classifierParams)
        elif self.classifierName =="logReg":
            return self.logReg(trainFeatures, trainClassLabels, evalFeatures,self.classifierParams)
        elif self.classifierName =="randC":
            return self.randC(trainFeatures, trainClassLabels, evalFeatures,self.classifierParams)
        ###TODO: add another classifiers here


    def randC(self, trainFeatures, trainClassLabels, evalFeatures, params = "default"):

        ind_random = np.random.randint(trainClassLabels.shape[0], size=(evalFeatures.shape[0],))
        evalClassLabels = trainClassLabels[ind_random]
        return evalClassLabels

    def mYkNN(self, trainFeatures, trainClassLabels, evalFeatures,params = "default"):
        """
        This is a varied version of kNN classifiers in sklearn used for DTW + KNN type classification, basically when the input features are distances already instance of features
        """
        #you have to estimate which is the test vector here, because based on that you have to estimate which column to take for kNN evaluation
        ind_eval = np.where(evalFeatures==0)[1]
        mYkNNParam = self.defaultClassifierSettings("kNN")
        if isinstance(params,dict):
            for key in params.keys():
                mYkNNParam[key]=params[key]

        mYkNNHandle = KNeighborsClassifier(**mYkNNParam)
        mYkNNHandle.fit(trainFeatures[:,ind_eval], trainClassLabels)
        evalClassLabels = mYkNNHandle.predict(evalFeatures[:,ind_eval])

        """
        THis is heuristic way of kNN which was kept to just crosschck the results
        K =5
        ind = np.argsort(trainFeatures[:,fold_ind])
        labels_kNN = trainClassLabels[ind[:K]]
        evalClassLabels1 = np.array([int(np.round((np.sum(labels_kNN)/5.0)))])

        if evalClassLabels[0]!=evalClassLabels1[0]:
            print "Here is the problem"
        """

        return evalClassLabels

    def nbMulti(self, trainFeatures, trainClassLabels, evalFeatures,params = "default"):

        nbMultiParam = self.defaultClassifierSettings("nbMulti")
        if isinstance(params,dict):
            for key in params.keys():
                nbMultiParam[key]=params[key]
        nbMultiHandle = MultinomialNB(**nbMultiParam)
        nbMultiHandle.fit(trainFeatures, trainClassLabels)
        evalClassLabels = nbMultiHandle.predict(evalFeatures)

        return evalClassLabels

    def logReg(self, trainFeatures, trainClassLabels, evalFeatures,params = "default"):

        logRegParam = self.defaultClassifierSettings("logReg")
        if isinstance(params,dict):
            for key in params.keys():
                logRegParam[key]=params[key]
        logRegHandle = LogisticRegression(**logRegParam)
        logRegHandle.fit(trainFeatures, trainClassLabels)
        evalClassLabels = logRegHandle.predict(evalFeatures)

        return evalClassLabels

    def kNN(self, trainFeatures, trainClassLabels, evalFeatures,params = "default"):

        kNNParam = self.defaultClassifierSettings("kNN")
        if isinstance(params,dict):
            for key in params.keys():
                kNNParam[key]=params[key]

        kNNHandle = KNeighborsClassifier(**kNNParam)
        kNNHandle.fit(trainFeatures, trainClassLabels)
        evalClassLabels = kNNHandle.predict(evalFeatures)

        return evalClassLabels

    def svm(self, trainFeatures, trainClassLabels, evalFeatures, params = "default"):
        """
        trains a svm classifier on trainFeatures using trainClassLabels as ground truths and then predict the class for features in evalFeatures
        Returns classLabels of the evalFeatures
        params is a dictionary to specify the parameters used in svm classifier in sklearn library. if its set to "default" (string) it will use default parameters.
        """

        svmParams = self.defaultClassifierSettings("svm")
        if isinstance(params,dict):
            for key in params.keys():
                svmParams[key]=params[key]

        svmHandle = SVC(**svmParams)
        #svmHandle = SVC(C=svmParams['C'], cache_size=svmParams['cache_size'], class_weight=svmParams['class_weight'], coef0=svmParams['coef0'], degree=svmParams['degree'], gamma=svmParams['gamma'], kernel=svmParams['kernel'], max_iter=svmParams['max_iter'], probability=svmParams['probability'], shrinking=svmParams['shrinking'], tol=svmParams['tol'], verbose=svmParams['verbose'])
        svmHandle.fit(trainFeatures, trainClassLabels)
        evalClassLabels = svmHandle.predict(evalFeatures)

        return evalClassLabels


    def tree(self, trainFeatures, trainClassLabels, evalFeatures, params = "default"):
        """
        trains a svm classifier on trainFeatures using trainClassLabels as ground truths and then predict the class for features in evalFeatures
        Returns classLabels of the evalFeatures
        params is a dictionary to specify the parameters used in svm classifier in sklearn library. if its set to "default" (string) it will use default parameters.
        """

        treeParams = self.defaultClassifierSettings("tree")
        if isinstance(params,dict):
            for key in params.keys():
                treeParams[key]=params[key]

        treeHandle = tree.DecisionTreeClassifier(**treeParams)
        #treeHandle = tree.DecisionTreeClassifier(criterion=treeParams["criterion"], max_depth=treeParams["max_depth"], min_samples_split=treeParams["min_samples_split"], min_samples_leaf=treeParams["min_samples_leaf"], min_density=treeParams["min_density"], max_features=treeParams["max_features"], compute_importances=treeParams["compute_importances"], random_state=treeParams["random_state"])
        treeHandle.fit(trainFeatures, trainClassLabels)
        evalClassLabels = treeHandle.predict(evalFeatures)

        return evalClassLabels
    
    def exportTREEModel(self, trainFeatures, trainClassLabels, modelFile, normFile):
        """
        this function exports a model built on trainingset
        output is stored in a pickle
        """
        #normalization of features, extracting mean variance of every dimension and saving
        norm = []
        for ii in np.arange(trainFeatures.shape[1]):
            norm.append({})
            norm[ii]['mean'] = float(np.mean(trainFeatures[:,ii]))
            norm[ii]['var'] = float(np.var(trainFeatures[:,ii]))
            trainFeatures[:,ii] = trainFeatures[:,ii]-norm[ii]['mean']
            trainFeatures[:,ii] = trainFeatures[:,ii]/norm[ii]['var']
        
        
        treeParams = self.defaultClassifierSettings("tree")
        treeParams['criterion']='entropy'
        treeHandle = tree.DecisionTreeClassifier(**treeParams)
        treeHandle.fit(trainFeatures, trainClassLabels)
        joblib.dump(treeHandle, modelFile)
        
        fid = file(normFile,'w')
        yaml.dump(norm,fid)
        fid.close()
        
    def predicByModel(self, modelFile, inpFeatures):
        
        model = joblib.load(modelFile)
        predictClassLabels = model.predict(inpFeatures)
        
        return predictClassLabels
        

    def defaultClassifierSettings(self, classifier):

        if classifier == "svm":
            return {"C":1.0, "cache_size":200, "class_weight": None, "coef0":0.0, "degree":3, "gamma":0.0, "kernel":'rbf', "max_iter":-1, "probability": False, "shrinking": True, "tol":0.001, "verbose":False}

        if classifier == "tree":
            return {"criterion":'gini', "max_depth": None, "min_samples_split": 2, "min_samples_leaf": 1, "min_density": 0.1, "max_features": None, "compute_importances":False, "random_state": None}

        if classifier == "kNN":
            return {"n_neighbors":5, "weights":'uniform', "algorithm":'auto', "leaf_size":30, "p":2, "metric":'minkowski'}

        if classifier == "nbMulti":
            return {"alpha":1.0, "fit_prior":True, "class_prior":None}

        if classifier == "logReg":
            return {"penalty":'l2', "dual":False, "tol":0.0001, "C":1.0, "fit_intercept":True, "intercept_scaling":1, "class_weight":None, "random_state":None}


    def setClassifierName(self, classifierName):
        self.classifierName = classifierName

    def setClassifierParams(self,classfierParams):
        self.classifierParams = classfierParams

    def requireNormFeatures(self,classifierName):

        if classifierName == "tree":
            return False
        if classifierName == "nbMulti":
            return False
        if classifierName == "svm":
            return True
        if classifierName == "kNN":
            return True
        if classifierName == "logReg":
            return True
        else:
            return True




class experimenter(classifiers):
    
    taskTitle = ""
    fNames = []
    cNames = []
    features = np.array([])
    classLabels = []
    classLabelsInt = []
    nExp = 0
    typeEval = ("kFoldCrossVal",10)
    
    def __init__(self):
        print "Initializing mlUtil class for using wrapper functions around core SkLearn library to perform experiments."
        #self.setExperimentParams()
        #classifiers.__init__(classifierName=self.classifierName, classifierParams=self.classifierParams)
        
    def readArffFile(self, filename):
        self.fNames = []
        f = open(filename,'r')        
        line = f.readline().strip()
        
        while line.lower() != "@data":
            linecontent = line.split()
            print linecontent
            if len(linecontent)>0 and linecontent[0] == "@relation":
                self.taskTitle = linecontent[1]
            
            if len(linecontent)>0 and linecontent[0] == "@attribute":
                
                if linecontent[2].strip() == "numeric":
                    self.fNames.append(linecontent[1])
                if linecontent[2].strip() != "numeric":
                    if linecontent[1] == "class":
                        classnames = ''.join(linecontent[2:])   # this is because sometimes class labels have white spaces
                        classnames = classnames.strip('{')
                        classnames = classnames.strip('}')
                        classnames = classnames.replace('\t','')
                        classnames = classnames.replace(' ','')
                        self.cNames = classnames.split(',')
                        self.cNames.remove('')  #this is because there are some extra ","s
                    else:
                        print "there is a problem reading %s file, features are not numeric" %filename
            line = f.readline().strip()
        
        contents = f.readlines()
        contents = [(j[:-1],j[-1]) for j in [ i.strip().split(",") for i in contents]]
        
        features = [f[0] for f in contents]
        classes = [c[1] for c in contents]
        
        self.features = np.array(features).astype(np.float)
        self.classLabels = np.array(classes)
        
        self.convertClassLabels2Int()
    
            
    def convertClassLabels2Int(self):
            self.classLabelsInt = copy.deepcopy(self.classLabels)

            for i, c in enumerate(self.cNames):
                ind = np.where(self.classLabelsInt == c)[0]
                self.classLabelsInt[ind] = i
            
            self.classLabelsInt= (self.classLabelsInt).astype(np.int)

    def setExperimentParams(self, nExp = 10, typeEval = ("kFoldCrossVal",10), nInstPerClass = -1, classifier = ('svm',"default")):
        self.nExp = nExp
        self.typeEval = typeEval
        self.nInstPerClass = nInstPerClass  #this means balance the dataset taking least number of instances present in a class

        self.setClassifierName(classifier[0])
        self.setClassifierParams(classifier[1])

        print "Number of experiments: %d \n"%self.nExp
        print "Type of evaluation: %s \n"%self.typeEval[0]
        print "Evaluation parameter for %s is: %d\n"%(self.typeEval[0], self.typeEval[1])
        print "Instances per class selected: %d\n" %self.nInstPerClass
        print "Classifier selected: %s\n"%self.classifierName


    def setFeaturesAndClassLabels(self, features, classLabels):

        self.features = features
        self.classLabels = classLabels
        self.cNames = list(set(classLabels))
        self.fNames = range(0,features.shape[1])
        self.convertClassLabels2Int()

    def genConfMTX(self, origClass, predClass):
        origClass = np.array(origClass)
        predClass = np.array(predClass)
        mtx = np.zeros([len(self.cNames),len(self.cNames)])

        for i in range(predClass.shape[0]):
            mtx[origClass[i],predClass[i]]+=1

        return mtx

    def normFeatures(self,array):

        return (array-np.mean(array))/np.sqrt(np.var(array))

    def runExperiment(self, features = -1, classLabels = -1, features2Use=-1):

        if not isinstance(features, int) and not isinstance(classLabels, int):
            self.setFeaturesAndClassLabels(features, classLabels)

        if self.features.shape[0] ==0:
            print "Please input features and class labels or provide path for arff file"
            return -1
        if len(self.classifierName) ==0:
            print "Before proceeding setExperimentParams()"
            return -1

        #Gathering all the indices corresponding to each class to further proceed for resampling dataset
        classMapp = []
        minInstPerClass = sys.float_info.max
        for i in range(0,len(self.cNames)):
            classMapp.append(np.where(self.classLabelsInt == i)[0])  #the indices of this classMapp array is the class names as the class name to Integer mapping is 0-N.. So classMap stores array of indices of members of that class
            print "Number of instances in class %s is %d \n"%(self.cNames[i],len(classMapp[i]))
            print type(classMapp[0])
            if(len(classMapp[i])<minInstPerClass):
                minInstPerClass = len(classMapp[i])
        


        ### Resampling of the dataset (if it should have equal number of features from both the classes (given a number N or minimum number of features of two classes)
        ### TODO: include sampling of instances of class which has few instances and repeat a random subset to increase the number of instances

        InstPerClass = np.zeros(len(self.cNames),dtype=int)
        if self.nInstPerClass == -1: #balancing the dataset with least number of instances present in any class
            InstPerClass[:] = minInstPerClass
        else:
            InstPerClass[:] = self.InstPerClass

        ### Variable to store confusion matrix of all experiments
        self.cMTXWhole = np.zeros([len(self.cNames),len(self.cNames)])
        self.cMTXExp = [[self.cMTXWhole] for x in range(self.nExp)]    #store confusion matrix of each iteration
        self.accuracy = [[] for x in range(self.nExp)]
        self.decArray = [[] for x in xrange(len(self.classLabelsInt))]  #stores all the decision for each feature vector in all the experiments


        #selection of the features
        ### Select all the user specified features
        if isinstance(features2Use,int):
            self.featuresSelected = copy.deepcopy(self.features)
        elif isinstance(features2Use, list):
            ind_features = []
            for feat in features2Use:
                ind_features.append(self.fNames.index(feat))
            self.featuresSelected = copy.deepcopy(self.features[:,ind_features])
        else:
            print "Features2Use should be either list of features or -1"

        # perform feature normalization and other procesing only when classifier is not mYkNN
        if not (self.classifierName == "mYkNN"):
            #normalization of the features
            for i in range(self.featuresSelected.shape[1]):
                self.featuresSelected[:,i] = self.normFeatures(self.featuresSelected[:,i])

            #making every feature non negative by giving DC offset
            if not self.requireNormFeatures(self.classifierName):
                for i in range(self.featuresSelected.shape[1]):
                    self.featuresSelected[:,i] = self.featuresSelected[:,i] + -1*np.min(self.featuresSelected[:,i])





        ### Experiment iterations start here
        for exp_cnt in range(0,self.nExp):

            expIndices = np.array([], dtype = int)  # this array will contain indices of instances which are to be used in one experiment
            for j in range(0,len(self.cNames)):     #reshuffling the order in every experiment before selecting the samples for each experiment
                np.random.shuffle(classMapp[j])
                expIndices = np.append(expIndices, classMapp[j][:InstPerClass[j]])

            ### Test and Evaluation part
            ### For each experiment perform test and evaluation in the following way.

            #shuffling indices of the chosen features for more randomness in selection of folds
            np.random.shuffle(expIndices)

            if (self.typeEval[1]==-1):
                n_folds = expIndices.shape[0]
            else:
                n_folds = self.typeEval[1]

            if self.typeEval[0] =="kFoldCrossVal":
                kfold = KFold(expIndices.shape[0], n_folds=n_folds, indices=True)

                for fold_ind,(train_ind, test_ind) in enumerate(kfold):
                    prediction = self.performTrainTest(self.featuresSelected[expIndices[train_ind],:],self.classLabelsInt[expIndices[train_ind]], self.featuresSelected[expIndices[test_ind],:])
                    resultPFold = np.where(prediction==self.classLabelsInt[expIndices[test_ind]])[0]

                    for k,ind in enumerate(test_ind):
                        self.decArray[expIndices[ind]].append(prediction[k])

                    self.accuracy[exp_cnt].append(float(resultPFold.shape[0])/float(test_ind.shape[0]))
                    cnfMTX = self.genConfMTX(self.classLabelsInt[expIndices[test_ind]],prediction)
                    self.cMTXExp[exp_cnt] = self.cMTXExp[exp_cnt] + cnfMTX
                    self.cMTXWhole = self.cMTXWhole + cnfMTX


            ### TODO: DO this for many type of evaluation methods

            
        self.overallAccuracy =[]
        for m in self.accuracy:
            self.overallAccuracy.append(np.mean(m))

        self.overallAccuracy = np.mean(self.overallAccuracy)

    def featureSelection(self,features2Use=-1):

        if isinstance(features2Use,int):
            self.featuresSelected = copy.deepcopy(self.features)
        elif isinstance(features2Use, list):
            ind_features = []
            for feat in features2Use:
                ind_features.append(self.fNames.index(feat))
            self.featuresSelected = copy.deepcopy(self.features[:,ind_features])
        else:
            print "Features2Use should be either list of features or -1"

    def normalizeFeatures(self):

        # perform feature normalization and other procesing only when classifier is not mYkNN
        if not (self.classifierName == "mYkNN"):
            #normalization of the features
            for i in range(self.featuresSelected.shape[1]):
                self.featuresSelected[:,i] = self.normFeatures(self.featuresSelected[:,i])

            #making every feature non negative by giving DC offset
            if not self.requireNormFeatures(self.classifierName):
                for i in range(self.featuresSelected.shape[1]):
                    self.featuresSelected[:,i] = self.featuresSelected[:,i] + -1*np.min(self.featuresSelected[:,i])




class advanceExperimenter(experimenter):

    def __init__(self, features=-1, labels=1, arffFile = -1):
        if isinstance(arffFile,str):
            self.readArffFile(arffFile)
        elif isinstance(features,np.ndarray):
            self.setFeaturesAndClassLabels(features, labels)

    def runCompoundExperiment(self, featureSets, classifierSets, experimentParams, resultDir, logFile):
        """
        This function performs compound experiments, for every featureset in FeatureSets and classifier in ClassifierSets running experiment with provided experimentParams. In the resultDir, logFile will contain a summary.
        also we store the result corresponding to each iteration and each fold whihc will be used for computing the significance test. This data is stored in a folder with the same name as logFile created inside resultDir
        """

        #create directories and file for dumping results and log of compund experiment
        out_dir = resultDir + '/' + (os.path.basename(logFile)).split('.')[0]+ "_outputs"
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.makedirs(out_dir)



        log = open(out_dir+"/"+logFile,'w')

        #dumping some information regarding the selected classifier set and feature set to mapp the indixes of experiments with the combination
        total_experiments = len(featureSets)*len(classifierSets)
        nFeatureSets = len(featureSets)
        nClassifier = len(classifierSets)

        log.write("This is the mapping of feature sets Vs their references in this log file:\n")
        for i, features in enumerate(featureSets):
            log.write(("featureSet"+str(i)+":"+ "\t%s"*len(features)+"\n")%tuple(features))

        log.write("\n")
        log.write("\n")
        log.write("\n")
        log.write("THis is the mapping of all the permutations of features with classifiers and experiment index:\n")
        log.write("\n")
        classifiers = [x[0] for x in classifierSets]
        log.write(("%s\t" + "\t%s"*len(classifierSets)+"\n")%tuple(["Ftr.\csf"]+classifiers))
        exp_cnt = 1
        for i,features in enumerate(featureSets):
            log.write("featureSet"+str(i)+"\t")
            for classifier in classifierSets:
                log.write("\t%d"%exp_cnt)
                exp_cnt+=1
            log.write("\n")

        log.write("\n")
        log.write("\n")
        log.write("\n")
        log.write("\n")

        self.setExperimentParams(**experimentParams)
        exp_cnt = 1
        log.write(("%s\t" + "\t%s"*len(classifierSets)+"\n")%tuple(["Ftr.\csf"]+classifiers))

        for i,features in enumerate(featureSets):
            log.write("featureSet"+str(i)+"\t")
            for classifier in classifierSets:
                self.setClassifierName(classifier[0])
                self.setClassifierParams(classifier[1])
                self.runExperiment(features2Use=features)
                fid = open(out_dir+"/"+str(exp_cnt)+".json",'w')
                json.dump(self.accuracy,fid, indent=4)
                fid.close()
                log.write("\t%0.3f"%self.overallAccuracy)
                exp_cnt+=1
            log.write("\n")



        log.close()

    def isStaticallyDifferent(self, result1, result2):

        result1= json.load(open(result1))
        result2= json.load(open(result2))

        result1 = np.reshape(result1, len(result1)*len(result1[0]))
        result2 = np.reshape(result2, len(result2)*len(result2[0]))
        U, p_val = mannwhitneyu(result1,result2)
        print U, p_val

        if p_val < .05:
            return True
        else:
            return False



