import re
#delete file before running script
#get data(ko file)

f = open('C:\Main\Cal\Research\kegg\module.txt')
data = f.read()
data = data.split('///')

#f = open('C:\Main\Cal\Research\kegg\module_sample.txt')
#data = f.read()



for module in data:
    lines = module.split('\n')

    ############################Get KOS########################
    for line in lines:
        if 'ORTHOLOGY' in line:
            start_key = lines.index(line)

        if 'REACTION' in line:
            stop_key = lines.index(line)
            break
        elif 'COMMENT' in line:
            stop_key = lines.index(line)
            break
        elif 'REFERENCE' in line:
            stop_key = lines.index(line)
            break
        
        


    kos = []
    for line in lines:
        if lines.index(line) >= start_key and lines.index(line) < stop_key:
            re_ko = re.compile('K[0-9,]+')
            #print line
            match = re_ko.search(line)
            if match:
                ko = match.group()
            if ko.find(',') != -1:
                ko_comma = ko.split(',')
                for ko_in in ko_comma:
                    kos.append(ko_in)
            else:
                kos.append(ko)

    #remove duplicates and bad kos
    kos = list(set(kos))
    for v in kos:
        re_num = re.compile('[0-9]')
        if re_num.search(v):
            print re_num.search(v).group
        else:
            print 'bad ko'
            kos.remove(v)

    #####################GET COMPOUNDS###########
    
    if module.find('COMPOUND') == -1:
        break
    
    for line in lines:
        if 'COMPOUND' in line:
            start_key = lines.index(line)
        
    cos = []

    for line in lines:
        if lines.index(line) >= start_key:
            re_co = re.compile('C[0-9]+')
            hit = re_co.search(line)
            if hit:
                co = hit.group()
                cos.append(co)

    #remove duplicates
    cos = list(set(cos))

    for a in cos:
        re_num = re.compile('[0-9]')
        if re_num.search(a):
            print re_num.search(a).group
        else:
            cos.remove(a)

    ####################GET ENTRY################

    re_en = re.compile ('M[0-9]+')

    for line in lines:
        if 'ENTRY' in line:
            entry_search = re_en.search(line)
            if entry_search:
                entry = entry_search.group()

    ##################INSERT INTO Comma Separated File#######

    for ko_entry in kos:
        entry_line = entry + ',' + ko_entry + "\n"
        filename = 'kos_module.csv'
        kf = open(filename,"a")
        kf.writelines(entry_line)

    for co_entry in cos:
        co_entry_line = entry + ',' + co_entry + "\n"
        filename = 'cos_module.csv'
        cf = open(filename,"a")
        cf.write(co_entry_line)


