class XY:
    
    def __init__(self):
        self.X = []
        self.Y = []

    def get_data(filename,alst,Complex):

        if Complex:
            with open(filename) as File:
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
                data = File.readline()
                for i in data.split(';'):
                    try:
                        alst.append(float(i))
                    except:
                        pass

    def get_X(self,filenameX,Complex=False):
        XY.get_data(filenameX,self.X,Complex)

    def get_Y(self,filenameY,Complex=False):
        XY.get_data(filenameY,self.Y,Complex)
    
