seg.extractSoloPercussionBATCHPROC('/media/Data/Datasets/MotifDiscovery_Dataset/CompMusic_Thodi/', '2s_MFCC_FLAT_CENT/2s_MFCC_FLAT_CENT_MODEL.pkl', '2s_MFCC_FLAT_CENT/2s_MFCC_FLAT_CENT_NORM.yaml',  0.0464399, 0.0266666666, 2, medianDur=60)

seg.exportTREEModel('2s_MFCC_FLAT_CENT/2s_MFCC_FLAT_CENT.arff', '2s_MFCC_FLAT_CENT/2s_MFCC_FLAT_CENT_MODEL.pkl', '2s_MFCC_FLAT_CENT/2s_MFCC_FLAT_CENT_NORM.yaml')

 seg.generateBinaryAggMFCCARFF('/media/Data/Datasets/percussionSegmentation/tani/', '/media/Data/Datasets/percussionSegmentation/nontani/', 'tani', 'non_tani', '2s_MFCC_FLAT_CENT.arff', 0.0464399, 0.0266666666, 2)