import numpy as np
import time
import os
import math
from matplotlib import pyplot as plt

def mean_cal(index,image):

   idx = index.copy()
   height = idx.shape[0]
   width = idx.shape[1]
   idx_matrix = np.zeros(shape=(height,width),dtype=np.uint32)
   label_index = []
   for h in range(height):
       for w in range(width):
           if idx[h, w] not in label_index:
               label_index.append(idx[h, w])
   idx_num = len(label_index)

   for h in range(height):
       for w in range(width):
           for i in range(idx_num):
               if (idx[h,w] == label_index[i]):
                   idx_matrix[h,w] = i

   adj_matrix = np.zeros(shape=(idx_num,idx_num),dtype= np.uint8)
   xh = [-1, -1, -1, 0, 0, 1, 1, 1]
   yh = [-1, 0, 1, -1, 1, -1, 0, 1]

   for h in range(height):
        for w in range(width):
            for i in range(8):
                if (h + xh[i] >= 0 and w + yh[i] >= 0 and h + xh[i] < height and w + yh[i] < width and
                    idx_matrix[h, w] != idx_matrix[h + xh[i], w + yh[i]]):
                    adj_matrix[int(idx_matrix[h, w])][int(idx_matrix[h + xh[i], w + yh[i]])] = 1

   img = image.copy()
   img_r = img[:,:,0]
   img_g = img[:,:,1]
   img_b = img[:,:,2]
   r_s = np.zeros(shape=(idx_num), dtype=np.float32)
   g_s = np.zeros(shape=(idx_num), dtype=np.float32)
   b_s = np.zeros(shape=(idx_num), dtype=np.float32)
   r_mean = np.zeros(shape=(idx_num),dtype=np.float32)
   g_mean = np.zeros(shape=(idx_num), dtype=np.float32)
   b_mean = np.zeros(shape=(idx_num), dtype=np.float32)
   num = np.zeros(shape=(idx_num), dtype=np.uint32)

   for i in range(idx_num):
       r = 0
       g = 0
       b = 0
       n = 0
       for h in range(height):
            for w in range(width):
               if idx_matrix[h,w] == i:
                   r += img_r[h,w]
                   g += img_g[h,w]
                   b += img_b[h,w]
                   n += 1
       r_s[i] = r
       g_s[i] = g
       b_s[i] = b
       num[i] = n
   for i in range(idx_num):
       r_mean[i] = r_s[i]/num[i]
       g_mean[i] = g_s[i]/num[i]
       b_mean[i] = b_s[i]/num[i]
   return adj_matrix, num, r_mean, g_mean, b_mean


def new_adj_process(adj_idx):

    num = adj_idx.shape[0]
    new_adj = np.zeros(shape=(num, num), dtype=np.uint8)

    for i in range(num):
        for j in range(num):
            if adj_idx[i, j] == 1 and adj_idx[j, i] == 1:
                new_adj[i, j] = 1
                new_adj[j, i] = 0
    return new_adj


def seg_show(image, index):
    m = index.shape[0]
    n = index.shape[1]
    image_new = image.copy()
    for h in range(1, m - 1):
        for w in range(1, n - 1):
            if (index[h,w] != index[h,  w+1] or index[h,w] != index[h, w-1] or
            index[h,w] != index[h+1, w] or index[h,w] != index[h-1,w]):
                image_new[h,w,:] = (255,0,0)
    return image_new


def reseg_idx(seg_result, index):
    num = len(seg_result)
    height = index.shape[0]
    width = index.shape[1]
    new_index = np.zeros(shape= (height, width), dtype= np.uint32)

    for i in range(num):
        new_id =  i
        ri = seg_result[i]
        r_num = len(ri)
        for j in range(r_num):
            id = ri[j]
            x_set = np.where(index == id)[0]
            y_set = np.where(index == id)[1]

            xn = len(x_set)
            for k in range(xn):
                dx = x_set[k]
                dy = y_set[k]

                new_index[dx, dy] = new_id

    return new_index


class graph_merging:
    def __init__(self,r_mean,g_mean,b_mean,num_matrix,adj_matrix,index,k):
        self.R_mean = r_mean
        self.G_mean = g_mean
        self.B_mean = b_mean
        self.R_mean = self.R_mean.astype(np.float16)
        self.G_mean = self.G_mean.astype(np.float16)
        self.B_mean = self.B_mean.astype(np.float16)
        self.count =  num_matrix
        self.adj_idx = adj_matrix

        self.idx = index

        self.num = self.count.shape[0]
        #self.time = tt
        self.k = k
        self.graph_processing()

    def graph_processing(self):
        self.new_adj = self.adj_idx


        self.result = []
        self.first_object()

        self.renew_object()




    def first_object(self):
        self.object=[]
        self.object_adj=[]

        for i in range(self.num):
            feature = [self.R_mean[i],self.B_mean[i],self.G_mean[i]]
            count = self.count[i]
            self.object.append((i,feature,count))
        #print((self.object[0])[1])
            adj_index =[]
            feature_adj = []
            count_adj = []
            for j in range(self.num):
                if self.new_adj[i,j] == 1:
                    adj_index.append(j)
                    feature_adj.append([self.R_mean[j],self.B_mean[j],self.G_mean[j]])
                    count_adj.append(self.count[j])
            self.object_adj.append((adj_index,feature_adj,count_adj))
        self.diff_first()
        #print(self.object_adj)
    def diff_first(self):
        self.diff_val= []
        for i in range(self.num):
            adj_num = len((self.object_adj[i])[0])
            #print(adj_num)
            feature_i = (self.object[i])[1]

            feature_adj = (self.object_adj[i])[1]
            diff = []

            for j in range(adj_num):
                a = self.distance(feature_i,feature_adj[j])
                diff.append(a)
            self.diff_val.append(diff)
        #print(self.diff_val)
        self.object_copy = self.object
        self.adj_copy = self.object_adj
        self.diff_copy = self.diff_val
        self.new_num = 0

    def distance(self, A, B):
        dif = [a - b for a, b in zip(A, B)]
        d = pow(dif[0], 2) + pow(dif[1], 2) + pow(dif[2], 2)

        return math.sqrt(d)

    def renew_object(self):
        self.Int = []
        self.one_list()
        ax = len(self.result)
        m = 0

    def one_list(self):
        adj_index = []
        diff_v = []
        num = len(self.object)
        for i in range(num):
            id = (self.object[i])[0]
            adj_num = len((self.object_adj[i])[0])
            idx = (self.object_adj[i])[0]
            dif = (self.diff_val[i])
            for j in range(adj_num):
                adj_index.append([id,idx[j]])
                diff_v.append(dif[j])

        a = np.sort(diff_v)
        b = np.argsort(diff_v)
        c = []
        for i in range(b.shape[0]):
            c.append(adj_index[b[i]])
        self.result.append(c[0])

        for i in range(1,b.shape[0]):
            r_num = len(self.result)
            for j in range(r_num):
                r_i = self.result[j]
                r_len = len(r_i)
                Int = []
                for r in r_i:
                    id = (self.object[r])[0]
                    id_adj = (self.object_adj[r])[0]
                    for rj in r_i:
                        if (rj in id_adj and rj != r):
                            r_a = id_adj.index(rj)
                            Int.append((self.diff_val[id])[r_a])
                if (Int == []):
                    Int2 = 0
                else:

                    Int2 = self.MInt(Int,r_len)
                self.Int.append(Int2)
            x_a = a[i]
            x_b = c[i]

            x1 = -1
            x2 = -1
            for j in range(r_num):
                r_i = self.result[j]
                r_len = len(r_i)
                for k in range(r_len):
                    if (r_i[k] == x_b[0]):
                        x1 = j

            for k in range(r_num):
                r_i = self.result[k]
                r_len = len(r_i)
                for k_2 in range(r_len):
                    if (r_i[k_2] == x_b[1]):
                        x2 = k

            if (x1 == -1 and x2 == -1):
                self.result.append(x_b)

            if (x1 != -1 and x2 == -1):
                x_s1 = self.Int[x1]
                if x_a <= x_s1:
                    self.result[x1].append(x_b[1])
                else:
                    self.result.append([x_b[1]])
            if (x1 == -1 and x2 != -1):
                    x_s2 = self.Int[x2]
                    if x_a <= x_s2:
                        self.result[x2].append(x_b[0])
                    else:
                        self.result.append([x_b[0]])

            if (x1 != -1 and x2 != -1):
                if (x1 == x2):
                    self.result = self.result
                elif (x1 != x2):
                    x_s1 = self.Int[x1]
                    x_s2 = self.Int[x2]
                    x_min = min(x_s1,x_s2)
                    #print(x_min)

                    if (x_a <= x_min):
                        r2_l = len(self.result[x2])
                        for rr in range(r2_l):
                            self.result[x1].append((self.result[x2])[rr])
                        self.result.remove(self.result[x2])
                    else:
                        self.result =self.result

            self.result = [x for x in self.result if x != []]

            #print(len(self.result))
            #print(len(self.Int))
            self.Int = []
        return self.result

    def MInt(self,Int,n):
        t = self.k / n
        b = [i + t for i in Int]
        x = max(b)

        return x

if __name__ == '__main__':
    main_folder = r"D:/pywork/large"
    file_str = ["ArcGIS", 'WDMI', 'google']

    for i in range(2,3):
        cur_str = file_str[i]

        start = time.time()

        file_folder = os.path.join(main_folder, cur_str)
        image_file_name = os.path.join(file_folder, cur_str + ".jpg")
        image = plt.imread(image_file_name)

        index_folder = os.path.join(file_folder, "Otsu_seg/mean/CMSuG")
        index_file = os.path.join(index_folder, "CMSuG_0_new.csv")

        index = np.loadtxt(open(index_file), delimiter=',', dtype = np.uint32, skiprows= 0)

        adj_matrix, num, r_mean, g_mean, b_mean = mean_cal(index,image)

        adj_new = new_adj_process(adj_matrix)

        k = 5

        seg_result = graph_merging(r_mean,g_mean,b_mean,num,adj_new,index,k)

        result = seg_result.result

        index_new = reseg_idx(result, index)

        tem = '%d' % k

        end = time.time()

        t = end - start

        save_folder = os.path.join(file_folder, "CMIS")
        save_csv = os.path.join(save_folder, "tau_" + str(tem) + ".csv")
        save_jpg = os.path.join(save_folder, "tau_" + str(tem) + ".jpg")
        image_new = seg_show(image, index_new)

        np.savetxt(save_csv, index_new, delimiter=',', fmt = '%d')
        plt.imsave(save_jpg, image_new)

        print(cur_str,"--", t)




