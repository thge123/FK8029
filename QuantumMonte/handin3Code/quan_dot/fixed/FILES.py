def get_data(filename,fileline=1,Complex=False):
    alst = []
    if Complex:
        with open(filename) as File:
            for i in range(fileline):
            	data = File.readline()
            for i in data.split(';'):
                z  = i[1:-1]
                z  = z.split(',')
                try:
                    alst.append((float(z[0]),float(z[1])))
                except:
                    pass
    else:
        with open(filename) as File:
            if fileline == 0:
                for i in File.readlines():
                    data = []
                    for j in i.split(';'):
                        try:
                            data.append(float(j))
                        except:
                            pass
                    alst.append(data)
                return alst
            else:
                for i in range(fileline):
                    data = File.readline()
                for i in data.split(';'):
                    try:
                        alst.append(float(i))
                    except:
                        pass
    return alst

def number2string(number):
    if number < 10:
        return '00'+str(number)
    elif 9 < number < 100:
        return '0'+str(number)
    elif 99 < number < 1000:
        return str(number)
    else:
        raise ValueError('Number invalid in number2string')

def string2number(string):
    return int(string[1:4])




