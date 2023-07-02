import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

import plantbox as pb
from rsml.rsml_data import RsmlData
import rsml.rsml_reader as rsml_reader
from estimate_params import *
import numpy as np
import os


class EstimateDataModel:
    """ 
    Model in the sense of model view controller (MVC), stores most that is presented in the view    
    """

    def __init__(self):
        self.rsmls = []  # list of rsml data to analyse, rsml_data.RsmlData
        self.estimates = [{}]  # list of dictionary of parameters per root
        self.plant = pb.Plant()  # to own the CPlantBox parameters
        self.parameters = []  # list of CPlantBox parameters of type , index = subType, max = 10
        self.pparameters = []  # list of CPlantBox Root System parameters
        self.orders = []  # list of root orders
        self.times = []
        self.base_root_indices = []  # tap roots and basal roots
        self.tap_root_indices = []
        self.basal_root_indices = []  # includes shoot bornes
        self.root_indices = []  # all roots
        self.folder_name = ""
        self.file_names = []  # for debugging

    def open_folder(self, folder_name):
        """ see RsmlData.open_rsml() in src/python_modules/rsml_data.py                
        """
        self.rsmls = []  # delete old data
        self.folder_name = folder_name
        for root, subdirs, files in os.walk(folder_name):
            for filename in files:
                if filename.endswith('.rsml'):
                    file_path = os.path.join(root, filename)
                    print()
                    print('file %s (full path: %s)\n' % (filename, file_path))
                    print()
                    self.file_names.append(filename)
                    file_data = RsmlData()
                    file_data.open_rsml(file_path)
                    if not "order" in file_data.properties:  # check if orders were created, if not create tag
                        file_data.properties["order"] = rsml_reader.get_root_orders(file_data.properties)
                    self.rsmls.append(file_data)
                    try:
                        if file_data.tagnames[1]:
                            self.times.append(file_data.max_ct)  #  check for creation_time tag
                            # self.times.append(11.)  # TESTING
                        else:
                            s = filename.split(".")
                            str = ''.join([n for n in s[0][-3:] if n.isdigit()])
                            self.times.append(int(str))
                    except:
                        print("filename", filename)
                        self.times.append(0)
        self.estimates = [None] * len(self.times)
        self.parameters = [pb.RootRandomParameter(self.plant) for _ in range(0, 10)]
        self.pparameters = [pb.SeedRandomParameter(self.plant)][0]
        self.create_length()  # add a length tag to properties
        self.initialize_roots_()  # find base roots

    def exists(self):
        """ true if a rsml folder was set """
        return len(self.rsmls) > 0

    def create_length(self):
        """ calculates root lengths, and adds it to properties if not already present)"""
        for k, _ in enumerate(self.rsmls):  # measurements
            if not "length" in self.rsmls[k].properties:
                self.rsmls[k].properties["length"] = []
                for i, p in enumerate(self.rsmls[k].polylines):
                    l = 0.
                    for j, n1 in enumerate(p[:-1]):
                        n2 = p[j + 1]
                        l += np.linalg.norm(np.array(n2) - np.array(n1))
                    self.rsmls[k].properties["length"].append(l)
        else:
            print("EstimateDataModel.create_length: 'length' tag is already available")

    def initialize_roots_(self):
        """ sets base root indices, tap root indices, and basal root indices.
        base root indices are the sum of taps and basals (treating shootborne and basals as equal for now) """
        # 1. find base root indices and all indices
        self.base_root_indices = []
        self.lat_root_indices = []
        self.root_indices = []
        for i, data in enumerate(self.rsmls):  # plants
            pp = data.properties["parent-poly"]
            self.base_root_indices.append([])
            self.lat_root_indices.append([])
            self.root_indices.append([])
            for j, _ in enumerate(data.polylines):
                self.root_indices[-1].append(j)
                if pp[j] == -1:  # add if base root
                    self.base_root_indices[-1].append(j)
                else:
                    self.lat_root_indices[-1].append(j)
                    # 2. find tap and basal indices
        n = len(self.times)
        self.tap_root_indices = [None] * n
        self.basal_root_indices = [None] * n
        for i in range(0, n):  # plants
            self.basal_root_indices[i] = []
            self.tap_root_indices[i] = []
            # print("measurement (file)", i)
            bri = np.array(self.base_root_indices[i], dtype = np.int64)
            lengths = np.array(self.rsmls[i].properties["length"])
            tap_index = np.argmax(lengths[bri])  # longest base root is tap root (not perfect)
            lo_ = np.argsort(lengths[bri])
            lo = lo_[::-1]
            self.base_root_indices[i] = [self.base_root_indices[i][j] for j in lo]
            self.tap_root_indices[i].append(self.base_root_indices[i][tap_index])
            # print("tap roots", self.tap_root_indices[-1], "tap index", tap_index, "root", self.base_root_indices[i][tap_index])
            self.basal_root_indices[i] = self.base_root_indices[i].copy()
            # print("base roots", self.base_root_indices[i], self.basal_root_indices[i])
            self.basal_root_indices[i].pop(tap_index)

    def create_params(self, apical_method, base_method, calibration_method):
        """
        """
        for i in range(0, len(self.times)):  # first make new estimates
            self.estimates[i] = {}

        # estimate zones per order
        order = 0
        indices = self.base_root_indices
        c = np.array([len(x) for x in indices])
        while np.sum(c) > 0:
            self.estimate_zones_(indices)
            # if base_method < 2:  # tap and basals are treated the same way (multiple dicots)
            #     self.aggregate_parameters_(indices, target_type = order)
            # else:
            #     pass  # TODO not implemented at many points in the code ####################################################################
            self.aggregate_parameters_(indices, target_type = order)
            order += 1
            indices = self.pick_order(order)  # update index set (TODO it must be always per order, but different target_types are possible for clustering and aggregation)
            c = np.array([len(x) for x in indices])

        # initial roots age
        for i, j_ in enumerate(self.base_root_indices):
            for j in j_:
                self.estimates[i][(j, "age")] = self.times[i]

        if base_method < 2:
            order = 0  # age must be estimated for intitial indices (!)
            indices = self.base_root_indices
            c = np.array([len(x) for x in indices])
        else:  # base_method >= 2
            order = 0
            indices = self.tap_root_indices  # tap root first
            self.fit_root_length_(indices, base_method, target_type = order)  # 1
            self.add_r_(indices, order)  # 2
            if apical_method == 0:  # delay based
                self.add_delay_(indices, order)
            self.compute_age(indices, order, apical_method)  # age for the next iteration, aggregate_parameters_ must be called before
            # copy r into basal parameters
            target_type = 4
            self.parameters[target_type] = self.parameters[0].copy(self.plant)
            # fit basal production rate
            indices = self.basal_root_indices
            p = self.parameters[target_type]
            n = len(indices)  # number of measurements
            length_basals = [None] * n
            for i, j_ in enumerate(indices):
                length_basals[i] = []
                for j in j_:
                    length_basals[i].append(self.rsmls[i].properties["length"][j])

            length_basals_ = [l or [1e-14] for l in length_basals]
            res, f, ages = estimate_order0_rate(np.array(length_basals_), p.r, p.lmax, self.times)
            if res.x[0] < 0:
                print("ERROR: negative growth rate!!!")
            else:
                print("production rate", res.x[0])

            # compute seed params
            delayB_, firstB_, maxB_ = [], [], []
            # print('times', self.times)
            for i in range(0, len(ages)):
                delayB_.append(np.mean(np.diff(np.sort(ages[i]))))
                firstB_.append(self.times[i] - np.max(ages[i]))
                maxB_.append(len(ages[i]))

            # print('delayB_, firstB_, maxB_', delayB_, firstB_, maxB_)

            srp = self.pparameters
            if delayB_:
                srp.delayB = np.nanmean(delayB_)
                srp.delayBs = np.nanstd(delayB_)
            if firstB_:
                srp.firstB = np.nanmean(firstB_)
                srp.firstBs = np.nanstd(firstB_)
            if maxB_:
                srp.maxB = np.nanmean(maxB_)
                srp.maxBs = np.nanstd(maxB_)
            # print('checkallParams', srp.delayB, srp.delayBs,srp.firstB,srp.firstBs, srp.maxB, srp.maxBs)

            # res, f, ages = estimate_order0_rrate(np.array(length_basals), p.r, p.lmax, self.times)
            # print("production rate", res.x[0], "elongation rate", res.x[1])

            # add basal ages accordingly
            # print('indices', indices)
            for i, j_ in enumerate(indices):
                for j in range(0, len(j_)):
                    # print('i,j', i,j)
                    # print('ages', ages)
                    self.estimates[i][(j_[j], "age")] = ages[i][j]
            self.fit_root_length_(indices, base_method, target_type = target_type)  # 1
            self.add_r_(indices, target_type)  # 2
            if apical_method == 0:  # delay based
                self.add_delay_(indices, target_type)
            self.compute_age(indices, target_type, apical_method)
            order = 1
            indices = self.pick_order(order)  # update index set (TODO it must be always per order, but different target_types are possible for clustering and aggregation)
            c = np.array([len(x) for x in indices])

        while np.sum(c) > 0:  # the rest (higher orders...)
            # 1. fit (r, lmax), based on age, into target_type = order
            self.fit_root_length_(indices, base_method, target_type = order)
            # 2. add individual r, for less noise
            self.add_r_(indices, order)
            if apical_method == 0:  # delay based
                self.add_delay_(indices, order)
            # 3. calculate next iteration ages
            self.compute_age(indices, order, apical_method)  # age for the next iteration, aggregate_parameters_ must be called before

            order += 1
            indices = self.pick_order(order)  # update index set (TODO it must be always per order, but different target_types are possible for clustering and aggregation)
            c = np.array([len(x) for x in indices])
            # print("new length", np.sum(c), "at order", order, "indices", c)  # TODO change while criteria [[],[],[]] will pass

    def estimate_zones_(self, indices):
        """ creates lb, ln, la per root (if possible)
        todo theta, a
        see aggregate_parameters, to obtain mean and sd 
        """
        for i, j_ in enumerate(indices):
            coordina = self.rsmls[i].polylines[0]
            for m in range(1, len(self.rsmls[i].polylines)):
                coordina = np.concatenate((coordina, self.rsmls[i].polylines[m]))
            for j in j_:
                kids = self.pick_kids(i, j)  # print("kids", kids[:, 0])
                # print("kids", len(kids[:, 0]))
                # print(kids.shape)
                if kids.shape[0] > 0:
                    ii = [self.rsmls[i].properties["parent-node"][kids[m, 0]] for m in range (0, kids.shape[0])]
                    base_polyline = self.rsmls[i].polylines[j]
                    ii.sort()
                    n0 = 0
                    nl0 = ii[0]
                    lb = polyline_length(n0, nl0, base_polyline)  # length of basal zone
                    ne = ii[-1]
                    nle = len(base_polyline) - 1
                    la = polyline_length(ne, nle, base_polyline)  # length of apical zone
                    ln_ = []
                    for k in range(0, len(ii) - 1):  # laterals
                        ln_.append(polyline_length(ii[k], ii[k + 1], base_polyline))  # internodal distances
                    # print("found values for ", i, j)
                    self.estimates[i][(j, "lb")] = lb
                    self.estimates[i][(j, "ln")] = ln_
                    self.estimates[i][(j, "la")] = la

                # set radius 'a'
                if "diameter" in self.rsmls[i].properties:
                    a_ = np.divide([self.rsmls[i].properties["diameter"][j]][0], 2)
                elif "diameter" in self.rsmls[i].functions:
                    a_ = np.mean(np.divide([self.rsmls[i].functions["diameter"][j]][0], 2))  # a is the mean of all roots of a specific order, all these roots are weighted equally, no matter how long they are
                else:
                    a_ = np.mean([self.rsmls[i].functions["radius"][j]][0])

                self.estimates[i][(j, "a")] = a_ / float(self.rsmls[i].metadata.resolution)

                # set branching angle 'theta'
                ii = [self.rsmls[i].properties["parent-node"]][0]
                if (int(ii[j]) == -1) or (int(ii[j]) == 1):  # basal root angle
                    if len(self.rsmls[i].polylines[j])>1: 
                        p = np.asarray(self.rsmls[i].polylines[j][0])
                        p2 = np.asarray(self.rsmls[i].polylines[j][1]) #take the node after since there is no previous node
                        v1 = [0, 0, -1]
                        v2 = p2 - p
                        if np.linalg.norm(v2)<1:
                            dummy = 2
                            while np.linalg.norm(v2)<1:
                                if len(self.rsmls[i].polylines[j])<=dummy:
                                    break
                                p2 = np.asarray(self.rsmls[i].polylines[j][dummy])
                                v2 = p2 - p
                                dummy = dummy+1
                    
                        v1 = v1 / np.linalg.norm(v1)
                        v2 = v2 / np.linalg.norm(v2)
                        #theta_ = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
                        angle = np.arccos(np.dot(v1, v2))/math.pi*180
                        #print('basal vectors', p, p2, angle)
                    else:
                        angle = 0   
                else:
                    p = coordina[(int(ii[j]))]
                    p0 = coordina[(int(ii[j] - 1))]
                    p2 = np.asarray(self.rsmls[i].polylines[j][0])
                    # print('p0, p, p2',p0, p, p2)

                    # growth angle should be measured between segments >1cm
                    v1 = p0 - p
                    if np.linalg.norm(v1) < 1:
                        dummy = 2
                        while np.linalg.norm(v1) < 1:
                            if ii[j] < dummy:
                                break
                            p0 = coordina[(int(ii[j] - dummy))]
                            v1 = p0 - p
                            dummy = dummy + 1

                    v2 = p2 - p
                    if np.linalg.norm(v2) < 1:
                        dummy = 1
                        while np.linalg.norm(v2) < 1:
                            if len(self.rsmls[i].polylines[j]) <= dummy:
                                break
                            p2 = np.asarray(self.rsmls[i].polylines[j][dummy])
                            v2 = p2 - p
                            dummy = dummy + 1

                    v1 = v1 / np.linalg.norm(v1)
                    v2 = v2 / np.linalg.norm(v2)
                    # theta_ = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
                    angle = np.arccos(np.dot(v1, v2)) / np.pi * 180

                if angle > 90:
                    angle = 180 - angle
                self.estimates[i][(j, "theta")] = angle

                # find maximum lateral root order
                if "order" in self.rsmls[i].properties:
                    order = [self.rsmls[i].properties["order"][j]][0]
                elif "order" in self.rsmls[i].functions:
                    order = np.max([self.rsmls[i].functions["order"][j]][0])
                elif "subType" in self.rsmls[i].properties:
                    order = [self.rsmls[i].properties["subType"][j]][0] - 1
                elif "subType" in self.rsmls[i].functions:
                    order = np.max([self.rsmls[i].functions["subType"][j]][0])
                else:
                    order = int(3)
                order = min(order, 3)
                self.estimates[i][(j, "order")] = order

    def aggregate_parameters_(self, indices, target_type):
        """ aggregates the individual root parameters into mean and sd, 
        and adds it to the parameters (of type list of RootRandomParameters) at index target_type  
        """
        la_, lb_, ln_, delay_, a_, theta_, order_ = [], [], [], [], [], [], []
        for i, j_ in enumerate(indices):
            for j in j_:
                # print('i,j', i,j)
                if (j, "la") in self.estimates[i]:
                    la_.append(self.estimates[i][(j, "la")])
                if (j, "delay") in self.estimates[i]:
                    delay_.append(self.estimates[i][(j, "delay")])
                if (j, "lb") in self.estimates[i]:
                    lb_.append(self.estimates[i][(j, "lb")])
                if (j, "ln") in self.estimates[i]:
                    ln_.extend(self.estimates[i][(j, "ln")])
                if (j, "a") in self.estimates[i]:
                    a_.append(self.estimates[i][(j, "a")])
                if (j, "theta") in self.estimates[i]:
                    theta_.append(self.estimates[i][(j, "theta")])
                if (j, "order") in self.estimates[i]:
                    order_.append(self.estimates[i][(j, "order")])

        p = self.parameters[target_type]
        if lb_:
            p.lb = np.mean(lb_)
            p.lbs = np.std(lb_)
        if ln_:
            p.ln = np.nanmean(ln_)
            p.lns = np.nanstd(ln_)
        if la_:
            # p.la = np.nanmean(la_)
            if np.isnan(p.ln):
                p.la = np.nanmean(la_)
            else:
                p.la = np.nanmean(la_) - p.ln / 2
            p.las = np.nanstd(la_)
        if delay_:
            p.ldelay = np.nanmean(delay_)
            p.ldelays = np.nanstd(delay_)
        p.a = np.nanmean(a_)
        p.a_s = np.nanstd(a_)
        p.theta = np.nanmean(theta_)
        p.thetas = np.nanstd(theta_)

        self.orders.append(np.max(order_))

    def compute_age(self, indices, target_type, apical_method):
        """ assumes the target_type has already a given age, and a fitted r, lmax;  
            calculates the age estimate of lateral roots            
        """
        p = self.parameters[target_type]
        for i, j_ in enumerate(indices):
            measurement_time = self.times[i]
            ct_name = self.rsmls[i].tagnames[1]
            ct_name = None  # TESTING
            if ct_name:
                for j in j_:
                    kids = self.pick_kids(i, j)
                    if len(kids) > 0:  # add estimated age
                        for k in kids[:, 0]:  # kids
                            if ct_name in self.rsmls[i].functions:
                                age = measurement_time - self.rsmls[i].functions[ct_name][k][0]
                            else:
                                age = measurement_time - self.rsmls[i].properties[ct_name]
                            self.estimates[i][(k, "age")] = age
            else:
                for j in j_:
                    kids = self.pick_kids(i, j)
                    if len(kids) > 0:  # add estimated age
                        parent_ct = measurement_time - self.estimates[i][(j, "age")]
                        for k in kids[:, 0]:  # kids

                            r = self.estimates[i][(j, "r")]  # of parent root
                            # r = self.parameters[target_type].r  # mean of parent root (DON'T)
                            il = self.rsmls[i].properties["parent-node"][k]  #
                            bl = polyline_length(0, il, self.rsmls[i].polylines[j])  # latearal base length

                            if apical_method == 0:
                                age = negexp_age(max(bl, 0.), r, p.lmax)
                                delay = self.estimates[i][(j, "delay")]
                                root_age = measurement_time - (parent_ct + age + delay)
                            else:
                                # ln2 = np.sum(self.estimates[i][(j, "ln")]) / 2
                                # if np.isnan(ln2):
                                #     ln2 = self.parameters[target_type].ln / 2
                                # la = self.estimates[i][(j, "la")] - self.parameters[target_type].ln / 2  # of parent root (MAYBE with individual ln/2)
                                la = self.parameters[target_type].la  # - self.parameters[target_type].ln / 2  # mean of parent root
                                bl_la = min(bl + la, p.lmax * 0.99)  # lmax of parent!!!
                                age = negexp_age(max(bl_la, 0.), r, p.lmax)  # time the lateral emerged
                                root_age = max(measurement_time - (parent_ct + age), 0.)
                                if root_age < 0:
                                    print(i, j, k)
                                    print(target_type, "measurement_time", measurement_time, "parent_ct", parent_ct, "age({:g},{:g})={:g}".format(bl + la, r, age),
                                          "root_age", root_age, "lmax", p.lmax, "lateral length", self.rsmls[i].properties["length"][k])
                                # root_age = max(root_age, 0.)
                                # root_age = min(root_age, measurement_time)

                            self.estimates[i][(k, "age")] = root_age  # tap root et =  0

    def fit_root_length_(self, indices, base_method, target_type):
        """ adds an estimate using base_method of r and lmax into target_type 
        """
        length_, age_ = [], []
        for i, j_ in enumerate(indices):
            for j in j_:
                age = self.estimates[i][(j, "age")]
                l = self.rsmls[i].properties["length"][j]
                age_.append(age)  # for fitting
                length_.append(l)  # for fitting
        length_ = np.array(length_)
        age_ = np.array(age_)

        if base_method == 0 or base_method == 2:
            r, k, res = fit_taproot_rk(length_, age_)
            print("order", target_type, "r", r, "k", k, "res", res)
        elif base_method == 1 or base_method == 3:
            k = self.parameters[target_type].lmax
            r, res = fit_taproot_r(length_, age_, k)
            print("order", target_type, "r", r, "k", k, "res", res)

        self.parameters[target_type].r = r
        self.parameters[target_type].lmax = k

    def add_r_(self, indices, target_type):
        """ adds r as drawn from the distribution, (r, lmax) must be fitted first """
        p = self.parameters[target_type]
        for i, j_ in enumerate(indices):
            for j in j_:
                l = self.rsmls[i].properties["length"][j]
                t = self.estimates[i][(j, "age")]
                l = min(l, p.lmax * 0.99)  #
                r = negexp_rate(l, p.lmax, t)
                self.estimates[i][(j, "r")] = r

    def add_delay_(self, indices, target_type):
        """ adds the apical delay based on measured la, r, and lmax """
        p = self.parameters[target_type]
        for i, j_ in enumerate(indices):
            for j in j_:
                kids = self.pick_kids(i, j)  # print("kids", kids[:, 0])
                if kids.shape[0] > 0:
                    r = self.estimates[i][(j, "r")]  # of parent root
                    ii = [self.rsmls[i].properties["parent-node"][kids[m, 0]] for m in range (0, kids.shape[0])]
                    base_polyline = self.rsmls[i].polylines[j]
                    ii.sort()
                    n0 = 0
                    nl0 = ii[0]
                    lb = polyline_length(n0, nl0, base_polyline)  # length of basal zone
                    ne = ii[-1]
                    nle = len(base_polyline) - 1
                    l0 = polyline_length(0, ne, base_polyline)  # length of apical zone
                    l1 = polyline_length(0, nle, base_polyline)  # length of apical zone
                    age0 = negexp_age(l0, r, p.lmax)
                    age1 = negexp_age(l1, r, p.lmax)
                    self.estimates[i][(j, "delay")] = age1 - age0

    def pick_order(self, root_order, order_tag = "order"):
        """
        """
        indices = []
        for i, data in enumerate(self.rsmls):  # plants
            indices.append([])
            for j, _ in enumerate(data.polylines):
                o = data.properties[order_tag][j]
                if o == root_order:  # add if of right topological order
                    indices[-1].append(j)
        return indices

    def pick_kids(self, measurement_index, base_root_index):
        """ returns the indices of all roots that are directly connected with 
        root @param index base_root in measurement @param measurement_index
        """
        ppi = np.array(self.rsmls[measurement_index].properties["parent-poly"])
        # print(base_root_index)
        # print("ppi", ppi)
        kids_indices = np.array(np.argwhere(ppi == base_root_index * np.ones(ppi.shape)), dtype = np.int64)
        # print(kids_indices[:, 0])
        # if len(kids_indices[:, 0]):
        #     input()
        return kids_indices

    def write_parameters(self, filename):
        """
        """
        for i in range(0, np.max(self.orders) + 1):  # only writes the existing number of root orders
            p = self.parameters[i]
            p.name = "root order {:g}".format(i)
            p.subType = i
            p.organType = 2
            self.plant.setOrganRandomParameter(p)
        # check if fibrous rootsystem (TODO: should be done with base_method maybe...)
        brc = []
        for elem in self.base_root_indices:
            brc.append(len(elem))
        if np.mean(brc) > 1:  # if fibrous root system
            srp = self.pparameters
            self.plant.setOrganRandomParameter(srp)
        self.plant.writeParameters(filename)

