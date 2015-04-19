def find_res(file):
    en_site = {}
    
    with open(file,'r') as thing:
        for line in thing:
            Line = line.replace('\n','')
            for i in range(len(Line)-2):
                 if Line[i] + Line[i+1] + Line[i+2] == '<1>':
                     enzyme = Line.replace('<','').replace('1','').replace('>','')
                 if Line [i] + Line[i+1] + Line[i+2] == '<3>':
    
                     if Line [i+3] + Line[i+4] + Line[i+5] + Line[i+6] + Line[i+7]== 'Mycob':
                         species = Line.replace('<','').replace('3','').replace('>','')
                         allclear = True
                        
                     else:
                         enzyme = ''
                         allclear = False
                 
                 if Line [i] + Line[i+1] + Line[i+2] == '<5>' and allclear:
                     
                     site = Line.replace('<','').replace('5','').replace('>','')
                     
                     en_site[species] = enzyme + '   ' + site
                     
    with open('restriction_enzymes_and_sites.txt','w') as f:
        for key in en_site:
            f.write(key + '     ' + en_site[key] + '\n\n')
            
            
    return en_site
                    
                