import os
import pandas as pd                                                                                                                                                                                     
import numpy as np                     
                                         
beataml = pd.read_csv(str(Path(os.environ.get("BEATAML_DIR", "..")) / "beataml_waves1to4_norm_exp_dbgap.txt", sep="\t", index_col=0, nrows=100).select_dtypes(include='number')                            
eisfeld = pd.read_csv(str(Path(os.environ.get("ARHG_BASE_DIR", "..")) / "Ras Signature" / "eisfeld_expression.tsv", sep="\t", index_col=0, nrows=100).select_dtypes(include='number')
                                                                                                                                                                                                          
print("BeatAML expression values:")                                                                                                                                                                     
print(f"  min={beataml.values.min():.2f}, max={beataml.values.max():.2f}, mean={beataml.values.mean():.2f}")                                                                                            
                                                                                                                                                                                                          
print("Eisfeld expression values:")    
print(f"  min={eisfeld.values.min():.2f}, max={eisfeld.values.max():.2f}, mean={eisfeld.values.mean():.2f}") 

from scipy.stats import zscore                                                                                                                                                                            
                                                                                                                                                                                                            
  # Z-score each gene (column) independently
eisfeld_zscored = pd.DataFrame(                                                                                                                                                                           
      zscore(eisfeld, axis=0),                                                                                                                                                                              
      index=eisfeld.index,
      columns=eisfeld.columns                                                                                                                                                                               
  ) 

