def get_data(filename):
    X = []
    with open(filename) as File:
        for i in File.readlines():
            X.append([float(j) for j in i.split(';')])
    return X
