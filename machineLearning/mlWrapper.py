import numpy as np
import copy, sys, os
import shutil
import math
import yaml
from sklearn.externals import joblib
from sklearn import svm , ensemble, linear_model, tree, neighbors, naive_bayes
#from sklearn.ensemble import RandomForestClassifier
#from sklearn.tree import DecisionTreeClassifier
#from sklearn.neighbors import KNeighborsClassifier
#from sklearn.naive_bayes import MultinomialNB
#from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing
import sklearn.cross_validation as crossVal
from collections import Counter
from scipy.stats import mannwhitneyu
import json


class classifiers():

    classifierName = ""
    classifierParams = {}
    classifierHandle = None
    skl_classifiers = {
                'svm': {'name': 'SVC', 'handle':svm.SVC, 'feat_norm':'z-norm', "norm_feat_req":True},
               'svm-lin': {'name': 'LinearSVC', 'handle':svm.LinearSVC, 'feat_norm':'z-norm', "norm_feat_req":True},
               'sgd': {'name': 'SGDClassifier', 'handle':linear_model.SGDClassifier, 'feat_norm':'z-norm', "norm_feat_req":True},
               'nb-multi': {'name': 'MultinomialNB', 'handle':naive_bayes.MultinomialNB, 'feat_norm': '', "norm_feat_req":False},
               'nb-gauss': {'name': 'GaussianNB', 'handle':naive_bayes.GaussianNB, 'feat_norm':'', "norm_feat_req":False},
               'nb-bern': {'name': 'BernoulliNB', 'handle':naive_bayes.BernoulliNB, 'feat_norm':'', "norm_feat_req":False},
               'tree': {'name': 'DecisionTreeClassifier', 'handle':tree.DecisionTreeClassifier, 'feat_norm':'', "norm_feat_req":False},
               'randForest': {'name': 'RandomForestClassifier', 'handle':ensemble.RandomForestClassifier, 'feat_norm':'', "norm_feat_req":False},
               'logReg': {'name': 'LogisticRegression', 'handle':linear_model.LogisticRegression, 'feat_norm':'z-norm', "norm_feat_req":True},
               'kNN': {'name': 'KNeighborsClassifier', 'handle':neighbors.KNeighborsClassifier, 'feat_norm':'z-norm', "norm_feat_req":True},
               'mYkNN': {'name': 'KNeighborsClassifier', 'handle':neighbors.KNeighborsClassifier, 'feat_norm':'z-norm', "norm_feat_req":True},
               }
    
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
        
        if not isinstance(classifierParams,dict):
            classifierParams = {}
        self.classifierParams = classifierParams
        
        if self.classifierName == "mYkNN":
            return self.mYkNN(trainFeatures, trainClassLabels, evalFeatures, self.classifierParams)
        elif self.classifierName == "randC":
            return self.randC(trainFeatures, trainClassLabels, evalFeatures, self.classifierParams)
        elif not self.skl_classifiers.has_key(self.classifierName):
            print "Please specify a valid classifier"
            return False
            
        classObj = self.skl_classifiers[self.classifierName]['handle'](**self.classifierParams)
        classObj.fit(trainFeatures, trainClassLabels)
        evalClassLabels = classObj.predict(evalFeatures)
        
        return evalClassLabels

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
        if not isinstance(params,dict):
            params = {}

        mYkNNHandle = self.skl_classifiers[self.classifierName]['handle'](**params)
        mYkNNHandle.fit(trainFeatures[:,ind_eval], trainClassLabels)
        evalClassLabels = mYkNNHandle.predict(evalFeatures[:,ind_eval])

        return evalClassLabels
    
    def exportClassifierModelHandle(self, trainFeatures, trainClassLabels, classifier, params):
        
        if not isinstance(params,dict):
            params = {}
        
        if not elf.skl_classifiers.has_key(classifier):
            print "Please provide a valid classifier"
            return False
        
        classHandle = self.skl_classifiers[classifier](**params)
        classHandle.fit(trainFeatures, trainClassLabels)
        
        return classHandle
    

    def setClassifierName(self, classifierName):
        self.classifierName = classifierName
    
    def setClassifierParams(self,classfierParams):
        self.classifierParams = classfierParams


class experimenter(classifiers):
    
    taskTitle = ""
    fNames = []
    cNames = []
    features = np.array([])
    classLabels = []
    classLabelsInt = []
    nExp = 0
    typeEval = ("kFoldCrossVal",10)
    ftr_scaler = {}
    
    def __init__(self):
        print "Initializing mlUtil class for using wrapper functions around core SkLearn library to perform experiments."
        self.ftr_scaler['z-norm'] = preprocessing.StandardScaler()
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
                        #self.cNames.remove('')  #this is because there are some extra ","s
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

    def setExperimentParams(self, nExp = 10, typeEval = ("kFoldCrossVal",10), nInstPerClass = -1, classifier = ('svm',"default"), balanceClasses=1):
        self.nExp = nExp
        self.typeEval = typeEval
        self.nInstPerClass = nInstPerClass  #this means balance the dataset taking least number of instances present in a class
        self.setClassifierName(classifier[0])
        self.setClassifierParams(classifier[1])
        self.balanceClasses = balanceClasses
        
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
    
    
    ############################ MAIN FUNCTION FOR RUNNING EXPERIMENTS ###########################
    
    def runExperiment(self, features = -1, classLabels = -1, features2Use=-1, verbose=1):
        """
        This is the main function to run experiments. There are two option to load data, either pass them as ndarray or use arff file (in which case use readArff function).
        
        Steps in this function:
        1) Handle Data IO (either through ndarray or through arff file)
        
        """
        
        #printing useful information 
        if verbose:
            print "Number of experiments: %d \n"%self.nExp
            print "Evaluation parameter for %s is: %d\n"%(self.typeEval[0], self.typeEval[1])
            print "Instances per class selected: %d\n" %self.nInstPerClass   
            print "Classifier selected: %s\n"%self.classifierName

        
        
        ### 1. Handle Data IO          (only if its provided in terms of ndarray). If you want to use arff file use readarff function before calling this function.
        if not isinstance(features, int) and not isinstance(classLabels, int):
            self.setFeaturesAndClassLabels(features, classLabels)
       
       
        ### Checking if there is not data
        if self.features.shape[0] ==0:
            print "Please input features and class labels or provide path for arff file"
            return -1
        
        
        ### Before trying to run experiment provide experiment parameters
        if len(self.classifierName) ==0:
            print "Before proceeding setExperimentParams()"
            return -1
        if not self.skl_classifiers.has_key(self.classifierName):
            print "please select a classifeir which is supported by this class"
            return -1
        
        if self.typeEval[0] =='leaveOneID':
            if isinstance(self.filterArray,int):
                print "In case of leaveOneID evaluation please provide filterArray containing ID mapping <use a function in this class to generate this array>"
                return -1
        ### TODO <Put additional checks here for other necessary parameters needed to run the experiment>

        ### Selection of the features (Select all the user specified features)
        if isinstance(features2Use,int):
            self.featuresSelected = copy.deepcopy(self.features)
        elif isinstance(features2Use, list):
            ind_features = []
            for feat in features2Use:
                ind_features.append(self.fNames.index(feat))
            self.featuresSelected = copy.deepcopy(self.features[:,ind_features])
        else:
            print "Features2Use should be either list of features or -1"
        
        # Soon after setting experimental params we can initialize variables to store the results/ other stats for each experiment
        ### Variable to store stats of all the experiments
        self.cMTXWhole = np.zeros([len(self.cNames),len(self.cNames)])
        self.cMTXExp = [copy.deepcopy(self.cMTXWhole) for x in range(self.nExp)]    #store confusion matrix of each iteration
        self.accuracy = [[] for x in range(self.nExp)]
        self.decArray = [[] for x in xrange(len(self.classLabelsInt))]  #stores all the decision for each feature vector in all the experiments
        

        
        ### FOLD GENERATION LOGIC is a crucial component of this functions there are various possibilities
        ## CASE 1) when filtering array is provided. Meaning no two token of the same id in filtering array can be in training and testing. Options here are:
        #               1) Leave one "id" cross validation (where one "id" is all the tokens of the same id in filtering array )
        #                       1) Balance the classes in training set or not. In this case there is no balance of classes in testing fold (since its determined by "id")
        #               2) Decide fold based on min # id count. 1 fold size = min(#"id"). IF YOU THINK DEEPLY SINCE THE TRAINING SET REMAINS SAME FOR EVREY FOLD THE CASE IS EXACTLY SIMILAR TO THE 1) in this category
        
        # CASE 2) When filtering array is not provided.
        #               1) K fold cross validation by having balanced classes in both testing and training set or not.
        
        
        ### Experiment iterations start here
        for exp_cnt in range(0,self.nExp):
            
            train_test_ind=[]

            # depending on the type of evaluation and its parameters (n_fold and balancing classes or not) determining training and testing indexes
            if self.typeEval[0] =='leaveOneID': #NOTE: that this is leave one ID (id is given to a bunch of features, eg. all the ones from a single file, do not confuse it with classes)
                
                expIndices = np.arange(self.featuresSelected.shape[0])  
                #counting unique ids
                idsArray = np.unique(self.filterArray)
                nInstPerClass = []
                for idfile in idsArray:
                    testInd = np.where(self.filterArray==idfile)[0] #reservind all the features for this specific file id for testing and the rest for training
                    restInd = np.setdiff1d(expIndices, testInd)
                    
                    # if the users have asked for balancing the dataset (equal number of instances per clases), this this evaluation mode, only the training set changes!
                    if self.balanceClasses:
                        classMapp = []
                        minInstPerClass = sys.float_info.max
                        for i in range(0,len(self.cNames)):
                            classind = np.where(self.classLabelsInt[restInd] == i)[0]
                            classMapp.append(restInd[classind])
                            if(len(classMapp[i])<minInstPerClass):
                                minInstPerClass = len(classMapp[i])
                        if self.nInstPerClass==-1:
                                nInstPerClass.append(minInstPerClass)
                        else:
                                nInstPerClass.append(self.nInstPerClass)
                        
                        trainInd = np.array([], dtype = int)  # this array will contain indices of instances which are to be used in one experiment
                        for j in range(0,len(self.cNames)):     #reshuffling the order in every experiment before selecting the samples for each experiment
                            np.random.shuffle(classMapp[j])
                            trainInd = np.append(trainInd, classMapp[j][:nInstPerClass[-1]])
                    else:
                        # if no balancing is required, then entire training set is used for training (no instances are removed because of balancing out)
                        trainInd = restInd
                    
                    train_test_ind.append((trainInd,testInd))

            elif  self.typeEval[0] == 'kFoldCrossVal' or self.typeEval[0] == 'kStratFoldCrossVal' or self.typeEval[0] == 'LeaveOneOut':
                # In the context of all these mode of evaluations, balancing classes means selecting equal number of feature vectors from every class.
                if self.balanceClasses:
                    classMapp = []
                    minInstPerClass = sys.float_info.max
                    for i in range(0,len(self.cNames)):
                        classMapp.append(np.where(self.classLabelsInt == i)[0])
                        if(len(classMapp[i])<minInstPerClass):
                            minInstPerClass = len(classMapp[i])
                    if self.nInstPerClass==-1:
                        self.nInstPerClass = minInstPerClass
                    
                    expIndices = np.array([], dtype = int)  # this array will contain indices of instances which are to be used in one experiment
                    for j in range(0,len(self.cNames)):     #reshuffling the order in every experiment before selecting the samples for each experiment
                        np.random.shuffle(classMapp[j])
                        expIndices = np.append(expIndices, classMapp[j][:self.nInstPerClass])
                else:
                    #if there is no balancing asked by the user then just take all the features
                    expIndices = np.arange(self.featuresSelected.shape[0])


                if self.typeEval[0] =='kStratFoldCrossVal':
                    if (self.typeEval[1]==-1):
                        n_folds = expIndices.shape[0]
                    else:
                        n_folds = self.typeEval[1]

                    np.random.shuffle(expIndices)
                    kfold = crossVal.StratifiedKFold(self.classLabelsInt[expIndices], n_folds=n_folds)
                    for train, test in kfold:
                        train_test_ind.append((train,test))
                
                elif self.typeEval[0] =='kFoldCrossVal':
                    if (self.typeEval[1]==-1):
                        n_folds = expIndices.shape[0]
                    else:
                        n_folds = self.typeEval[1]
                        
                    np.random.shuffle(expIndices)
                    kfold = crossVal.KFold(expIndices.shape[0], n_folds=n_folds)
                    for train, test in kfold:
                        train_test_ind.append((train,test))

                elif self.typeEval[0] =='LeaveOneOut':
                    kfold = crossVal.LeaveOneOut(expIndices.shape[0])
                    for train, test in kfold:
                        train_test_ind.append((train,test))           
            else:
                print "Please provide a valid evaluation type"

            for fold_ind,(train_ind, test_ind) in enumerate(train_test_ind):
                ftr_train = copy.deepcopy(self.featuresSelected[expIndices[train_ind],:])
                ftr_test = copy.deepcopy(self.featuresSelected[expIndices[test_ind],:])
                if self.skl_classifiers[self.classifierName]['norm_feat_req']:
                    ftr_train = self.ftr_scaler[self.skl_classifiers[self.classifierName]['feat_norm']].fit_transform(ftr_train)
                    ftr_test = self.ftr_scaler[self.skl_classifiers[self.classifierName]['feat_norm']].transform(ftr_test)
                prediction = self.performTrainTest(ftr_train,self.classLabelsInt[expIndices[train_ind]], ftr_test)
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
        
    def generateFilterArrayFromMappFile(self, mappFile):
        
        lines = open(mappFile,'r').readlines()
        nFeatArray = np.array([])
        for line in lines:
            nFeatures = line.split('\t')[-1].strip()
            nFeatArray = np.append(nFeatArray, nFeatures)
        
        nFeatArray = nFeatArray.astype(np.int)
        
        filterArray = np.ones(np.sum(nFeatArray))
        last_ind = 0
        for ii,increment in enumerate(nFeatArray):
            filterArray[last_ind:last_ind+increment]=ii
            last_ind = last_ind+increment
            
        
        self.filterArray = filterArray
        
        return filterArray

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

    def performFeatureNormalization(self):
        """
        This function performs feature normalization depending upon the kind of classifier that has to be used
        """
        normType = 'Null'
        
        if self.classifierName == "tree":
            normType = 'AllPositive'
        if self.classifierName == "nbMulti":
            normType = 'AllPositive'
        if self.classifierName == "svm":
            normType = 'ZeroMean1Var'
        if self.classifierName == "kNN":
            normType = 'ZeroMean1Var'
        if self.classifierName == "logReg":
            normType = 'ZeroMean1Var'
        if self.classifierName == "mYkNN":#this function is a hack function, see description before making judgements
            pass
        
        self.normFactors = []
        if normType == 'ZeroMean1Var':
            for i in range(self.featuresSelected.shape[1]):
                 mean = np.mean(self.featuresSelected[:,i])
                 std = np.std(self.featuresSelected[:,i])
                 self.normFactors.append((i,mean,std))
                 self.featuresSelected[:,i] = (self.featuresSelected[:,i]-mean)/std
        elif  normType == 'AllPositive':
            for i in range(self.featuresSelected.shape[1]):
                    min_val = np.min(self.featuresSelected[:,i])
                    self.normFactors.append((i,min_val,1))
                    self.featuresSelected[:,i] = self.featuresSelected[:,i] -min_val

        
    def exportModel(self, classifier, classifierParams, modelFile, normFile):
        """
        this function exports a model built on trainingset
        output is stored in a pickle
        """
        self.setClassifierName(classifier)
        self.setClassifierParams(classifierParams)
        self.featureSelection()
        self.performFeatureNormalization()
        
        classHandle = self.exportClassifierModelHandle(self.featuresSelected, self.classLabelsInt, self.classifierName, self.classifierParams)
        if isinstance(classHandle, int):
            if classHandle ==-1:
                print "Model couldn't be exported"
                return -1
            
        joblib.dump(classHandle, modelFile)
        
        norm = []
        for ii in np.arange(self.featuresSelected.shape[1]):
            norm.append({})
            norm[ii]['mean'] = self.normFactors[ii][1]
            norm[ii]['var'] = self.normFactors[ii][2]
            
        fid = file(normFile,'w')
        yaml.dump(norm,fid)
        fid.close()
        
        
    def predicByModel(self, modelFile, inpFeatures):
        
        model = joblib.load(modelFile)
        predictClassLabels = model.predict(inpFeatures)
        
        return predictClassLabels
                        




class advanceExperimenter(experimenter):

    def __init__(self, features=-1, labels=1, arffFile = -1, mappFile=-1):
        if isinstance(arffFile,str):
            self.readArffFile(arffFile)
        elif isinstance(features,np.ndarray):
            self.setFeaturesAndClassLabels(features, labels)
        
        if isinstance(mappFile, str):
            self.generateFilterArrayFromMappFile(mappFile)

    def runCompoundExperiment(self, featureSets, classifierSets, experimentParams, resultDir, logFile):
        """
        This function performs compound experiments, for every featureset in FeatureSets and classifier in ClassifierSets running experiment with provided experimentParams. In the resultDir, logFile will contain a summary.
        also we store the result corresponding to each iteration and each fold whihc will be used for computing the significance test. This data is stored in a folder with the same name as logFile created inside resultDir
        """

        #create directories and file for dumping results and log of compund experiment
        out_dir = resultDir
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



