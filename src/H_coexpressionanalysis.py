#before this make sure to'pip install PyWGCNA'
#note, I had to install rust and visual C++ prerequisites
#before installing PyWGCNA (not enough space tho)
print("I am running!")

import PyWGCNA
print("imported PyWGCNA")
import pandas as pd
import os
import colorsys
print("finished imports")

class coexpressionanalysis:

  def __init__(self):
    pass

  def WGCNA(self, geneExppaths, metadatapaths, variables, correlationtype, outputpath):
    #note the geneExpfilepaths, metadatapaths and variables are lists
    #of all the paths for the different files. One file per variable

    #First prepare the input to be in the right format:
    #rows/obs are samples/sale info and cols/var are genes/gene info
    
    #create a list to store the individual dataframes
    print("create list to store the individual dataframes")
    dfs = []
    for i in geneExppaths:
      #create database from csv file
      df = pd.read_csv(i, index_col=0)
      #in this dataframe, each row is a gene and each column is a sample 
      #append them to make one big dataframe
      dfs.append(df)
    #concatenate dataframes, but first check that indexes match
    if all(df.index.equals(dfs[0].index) for df in dfs):
      combined_df = pd.concat(dfs, axis=1)
    else:
      print('Sample row indices do not match, cannot concatenate')
      pass
    #transpose it to get it in the right format
    geneexp = combined_df.T
    #print('transformed df')     #troubleshooting
    #print(geneexp.head())   #troubleshooting
        
    ####This section is not necessary ***
    '''
    #making it a csv file to be used as an input for PyWGCNA
    combined_df.to_csv(path_or_buf=os.path.join(outputpath,'geneexpression.csv'))
    filepath = os.path.join(outputpath,'geneexpression.csv')
    '''
      
    #initialyze WGCNA object
    print("initialize WGCNA object")
    wgcna = PyWGCNA.WGCNA(geneExp=geneexp, outputPath=outputpath, save=True)
    print(f"Raw expresion data along with information:\n {wgcna.geneExpr}")
    wgcna.geneExpr.to_df().head(5)


    #preprocessing and module building
    print("preprocessing and module building")
    wgcna.runWGCNA()
    print(f"Processed expresion data along with information:\n {wgcna.datExpr}")
    print(f"Raw gene inforamtions:")
    print(wgcna.datExpr.var.head(5))

      
    #inserting metadata *****************Test this section*********
    #the variable order corresponds to the geneExppaths and metadatapaths order
    #giving each condition a color and separating by variables for later analysis
    #initiating empty list to store metadataframes
    metadfs = []
    #initializing color mapping dictionary
    color_mapping = {}
    print("concatenate metadfs")
    for i in range(len(metadatapaths)): 
      #create database from csv file
      metadf = pd.read_csv(metadatapaths[i], index_col=0)
      metadfs.append(metadf)
      unique_conditions = metadf['condition'].unique()
      
      '''print("setting metadata colours")
      for j, condition in enumerate(unique_conditions):
          # Generate color in HSV color space with equally spaced hues
          hue = (j * 1.0 / len(unique_conditions)) % 1.0
          # Convert HSV color to RGB
          rgb_color = colorsys.hsv_to_rgb(hue, 0.9, 0.9) #(hue, saturation, value)
          # Convert RGB values to hexadecimal string
          hex_color = '#{:02x}{:02x}{:02x}'.format(int(rgb_color[0] * 255),
                                                    int(rgb_color[1] * 255),
                                                    int(rgb_color[2] * 255))
          color_mapping[condition] = hex_color
          # Set metadata color using setMetadataColor()
          #wgcna.setMetadataColor(condition, {condition: hex_color})
      wgcna.setMetadataColor('condition', color_mapping) '''
        ### This will not work if there are multiple metadata paths, bc a new 
          ### colour mapping will be assigned to 'condition' on every iteration, 
          ### and only the last one will be kept
      '''for j in range(len(unique_conditions)):
        #creating dictionary of colors 
        color_mapping.update((unique_conditions[j],wgcna.metadataColors[j]))
      wgcna.setMetadataColor(variables[i], color_mapping)'''
    
    if all(metadf.index.equals(metadfs[0].index) for metadf in metadfs):
      combined_metadf = pd.concat(metadfs, axis=1)
    else:
      print('Metadata row indices do not match, cannot concatenate')
      pass
    
    print("update sample info in WGCNA object")
    wgcna.updateSampleInfo(sampleInfo=combined_metadf, sep=',')
    print('sample info')              #troubleshooting
    print(wgcna.geneExpr.obs.head(5)) #troubleshooting
    print(f"Raw gene inforamtions:")
    print(wgcna.datExpr.var.head(5))

    color_mapping = {}
    print("concatenate metadfs")
    for metadf_i in metadfs:
      #condition_header=metadf_i[0]
      unique_conditions = metadf_i['condition'].unique()
      print("setting metadata colours")
      for j, condition in enumerate(unique_conditions):
          # Generate color in HSV color space with equally spaced hues
          hue = (j * 1.0 / len(unique_conditions)) % 1.0
          # Convert HSV color to RGB
          rgb_color = colorsys.hsv_to_rgb(hue, 0.9, 0.9) #(hue, saturation, value)
          # Convert RGB values to hexadecimal string
          hex_color = '#{:02x}{:02x}{:02x}'.format(int(rgb_color[0] * 255),
                                                    int(rgb_color[1] * 255),
                                                    int(rgb_color[2] * 255))
          color_mapping[condition] = hex_color
          # Set metadata color using setMetadataColor()
          #wgcna.setMetadataColor(condition, {condition: hex_color})
      wgcna.setMetadataColor('condition', color_mapping) 
        ### This will not work if there are multiple metadata paths, bc a new 
          ### colour mapping will be assigned to 'condition' on every iteration, 
          ### and only the last one will be kept

    print("metadata colours")
    print(wgcna.metadataColors)
    print(wgcna.metadataColors.keys())

    '''
    #### we don't need this for our data but we may want to include it
    #updating gene info
    geneList = PyWGCNA.getGeneList(dataset='mmusculus_gene_ensembl',
                                  attributes=['ensembl_gene_id', 
                                              'external_gene_name', 
                                              'gene_biotype'],
                                  maps=['gene_id', 'gene_name', 'gene_biotype'])
    wgcna.updateGeneInfo(geneList)
    '''


    print("analyzing results")
    #analyzing results
    wgcna.analyseWGCNA()
    #save object as a .p file
    print("save WGCNA")
    wgcna.saveWGCNA()

    #PATHWAY CO-EXPRESSION

    #asking first for other geneIDs ***make sure the input is in the right format**
    pathwaygenes = input('Please input a list of geneIDs of genes you know are involved in the pathway of interest: ')
    # Find the intersection of the requested genes and the genes in the dataset
    valid_genes = set(pathwaygenes).intersection(set(wgcna.datExpr.var_names))
    #make a list of the modules that other pathway genes belong to
    pathwaymodules=[]
    for geneId in valid_genes:
        pathwaymodules.append(wgcna.datExpr.var.moduleColors[geneId])

    if len(pathwaymodules) == 1:
      pathwaymodules = pathwaymodules[0]

    pathwaymodules = list(set(pathwaymodules)) #keep only the unique modules

    wgcna.outputPath = "Results" ## THIS MIGHT NEED TO BE ADJUSTED 
      #(results from CoExpressionModulePlot are put into ./Results    I'm not sure about other commands)

    ###!!! I haven't confirmed that this part works properly yet
    #Get candidate genes for modules with high correlations to the variables of interest
    ##### if we are not assuming that all the data they give us is for variables of interest we need to add
    ##### an input that is a list of the variables of interest only
    #get the module trait correlations
    module_trait_correlation = wgcna.moduleTraitCor
    #this part will select the modules with the proper variable correlation from the pathway modules
      
    #initialize a properly-correlated genes dataframe
    correlatedgenesdf = pd.DataFrame(columns =['Variable', 'Module','Genes'])

    #look for modules that are correlated to the variables in the desired way
    for i in range(len(variables)):
      if correlationtype[i]=='positive':
        #iterate through the pathway modules
        for j in pathwaymodules:
          #if the correlation to the variable is positive
          if module_trait_correlation.loc[j,variables[i]]>0:
            correlatedgenes = wgcna.getGeneModule(j)
            correlatedgenesdf = correlatedgenesdf.append(
              {'Variable': variables[i],
                'Module': j,
                'Genes':', '.join(correlatedgenes)}, 
                ignore_index=True
              )
      elif correlationtype[i]=='negative':
        #iterate through the pathway modules
        for j in pathwaymodules:
          #if the correlation to the variable is positive
          if module_trait_correlation.loc[j,variables[i]]<0:
            correlatedgenes = wgcna.getGeneModule(j)
            correlatedgenesdf = correlatedgenesdf.append(
              {'Variable': variables[i],
                'Module': j,
                'Genes':', '.join(correlatedgenes)}, 
                ignore_index=True
              )
    ###!!!

    #now save data frame to csv file
    correlatedgenesdf.to_csv(os.path.join(outputpath,'candidate genes from pathway coexpression with experimental variable correlation'))

    ####visualization for the POSTER
    #for module visualization
    wgcna.CoexpressionModulePlot(modules=pathwaymodules, numGenes=10, numConnections=100, minTOM=0)
    #for correlation visualization
    wgcna.module_trait_relationships_heatmap(metadfs) # Requires a metadata entry but idk what the metadata is supposed to look like
    
  pass



#In the main file, execution would be something like this

directory = '/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/Results'
geneexpressionpaths = ['/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/References/var_1/vst_counts.csv']
metadatapaths = ['/Users/xcatherine/Desktop/BIEN470/bioinformatic-pipeline/References/var_1/metadata.csv']
variables = ['light']
correlationtype = ['negative']

#cea = coexpressionanalysis()
#cea.WGCNA(geneexpressionpaths, metadatapaths,variables, correlationtype, directory)




