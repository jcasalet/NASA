PROTOBUF_FILE=/Users/jcasalet/Desktop/FLUID/CRISP/data/WORKSPACE/EARTH/AGG/workspace/save/crisp_best_nlerm.pt
TRAINDATA_FILE=/Users/jcasalet/Desktop/FLUID/CRISP/data/col_0/train/data.csv



python get_results.py -pf $PROTOBUF_FILE -tf $TRAINDATA_FILE -nf 10 -mn ERM -pw False
