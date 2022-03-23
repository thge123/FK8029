class XY:
    
    def __init__(self):
        self.X = []
        self.Y = []

    def get_data(filename,alst):
        with open(filename) as File:
            data = File.readline()
            for i in data.split(';'):
                try:
                    alst.append(float(i))
                except:
                    pass

    def get_X(self,filenameX):
        XY.get_data(filenameX,self.X)

    def get_Y(self,filenameY):
        XY.get_data(filenameY,self.Y)
