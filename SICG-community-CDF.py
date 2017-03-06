# -*- coding:GBK -*-
from itertools import combinations 
import networkx as nx
import numpy as np
import igraph as igraph
#import scipy.special as special
#import scipy.optimize as opt
#import matplotlib.pyplot as plt       
import sys
import os
from math import *

nodes = []
edges = []
nodemap = {} #node id and nodename

writedir=sys.argv[2]
os.chdir(writedir)
def analysis_depen(filename):
    G = nx.DiGraph()
    findFile = open(filename,'r')
    package=''
    each_lines= findFile.readlines()
    for each_line in each_lines:
        if each_line.__contains__('->'):
	#说明是一条边的记录
            edge=each_line.split('->');
            w=edge[1][edge[1].index('[label=')+7:edge[1].index(',')]
            if(edge[1].__contains__('[')):
                edge[1]=edge[1][:edge[1].index('[')]
            edge[0]=nodemap[edge[0]]+'-'+edge[0]
            edge[1]=nodemap[edge[1]]+'-'+edge[1]
            G.add_edge(edge[0],edge[1],weight=float(w))
            #print(edge[0]+'到'+edge[1]) #调试语句，显示边连接关系
        else:
        #说明是一条点的记录
            if each_line.__contains__('[label='):
                key=each_line[each_line.index('\"')+1:each_line.index('[')-2]
                if key=='depNode_0':
                    node='$'
                else:
                    node=each_line[each_line.index('[label=')+8:each_line.rindex('(')]
                    if package!='':
                        node=package+'.'+node
                nodemap[key] = node
                G.add_node(node+'-'+key)
                '''
                由于存在函数名重载的情况，所以将方法名加上节点编号共同作为图中的节点名
                '''
            else:
                if each_line.__contains__('@'):
                    package=each_line[each_line.index(':')+1:each_line.rindex('\"')]
                else:
                    if each_line.__contains__('}'):
                        package=''#说明到了一个包的结尾处，此时将package置为空               
    findFile.close()
    print(str(len(nodemap)))
    G_new = G.copy() #define a new directed graph
    return G_new

def prob_density_function(xs,norm=False):
    distKeys = range(min(xs),max(xs)+1)
    pdf = dict([(k,0) for k in distKeys])
    for x in xs:
        pdf[x] += 1
    if norm:
        pdf.update([(k,nf) for (k,nf) in zip(pdf.keys(),[float(f)/sum(pdf.values()) for f in pdf.values()])])
    return pdf

def cum_density_function(xs,norm=False):
    pdf = prob_density_function(xs)
    cdf = dict([(k,cf) for (k,cf) in zip(pdf.keys(),np.cumsum(pdf.values()))])
    if norm:
        cdf.update([(k,ncf) for (k,ncf) in zip(cdf.keys(),[float(cf)/sum(pdf.values()) for cf in cdf.values()])])
    return cdf

def complement_cumulative_distribution(xs):
    """Return the complement cumulative distribution(CCD) of a list dist

    Returns the CCD given by
            P(X>x) = 1 - P(X<x)
    where P(X<x) is the cumulative distribution function given by
            P(X < x) = sum_{xi<x} p(xi)
    where p(xi) is the probaility density function calculated as
            p(xi) = xi/(sum_i xi)

    Parameters
    ----------
    dist : list of Numeric
           list of values representing the frequency of a value at
           the index occuring
    norm : Whether or not to normalize the values by the sum

    Returns
    -------
    ccdh : dict of floats, keyed by occurance
           A dict of the same length, as the cumulative complementary
           distribution function.
    """
    ccdf = cum_density_function(xs,norm=True)
    ccdf.update([(k,1-cf) for (k,cf) in zip(ccdf.keys(),ccdf.values())])
    return ccdf

def discrete_power_law_neg_log_likelihood(a,xmin,xs):
    return float(len(xs)*np.log(special.zeta(a,xmin)) + a*sum([np.log(x) for x in xs]))

def discrete_power_law_exponent_estimator(x,xmin=1,a0=2.5):
    """ Returns the exponent of the power law of the degree distribution
    or at least attempts to find one.

    If the degree distribution follows the power law

        p(x) \prop x^-a

    where 2 < a < 3 generally.[1]

        "... one can estimate a by direct numeric maximization of the
        likelihood function itself, or equivalently of its logarithm
        (which is usually simpler):

        L(a) = -n ln(zeta(a,xmin)) - a sum(ln(xi))"[1]

    in the discrete case.

    We actually minimize the the -log likelihood as this is the same

    While this function will actually return the exponent for the graph
    this doesn't necessarily mean the data has a power law
    distribution. Other statistical tests should be run to build evidence
    for that hypothesis.
    
    Parameters
    ----------
    x    : list of collected values
    xmin : minimum value for which to estimate power law value
    a0   : starting guess at exponent

    Returns
    ----------
    a : float
        The estimation of the scaling parameter

    [1] A. Clauset, C.R. Shalizi, and M.E.J. Newman,
        "Power-law distributions in empirical data"
        SIAM Review 51(4), 661-703 (2009). (arXiv:0706.1062)

    A little explaining about this function opt.fmin. It takes a
    function that represents the function we are trying to minimize,
    in this case the negative log likelihood of the discrete power
    law function. It takes an initial guess 2.5, which is where most
    power laws are anyway and args which tells it what arguments to
    pass to the function we are optimizing besides a, and finally we
    turn the display off with disp=0
    """
    a = opt.fmin(discrete_power_law_neg_log_likelihood,a0,args=(xmin,x),disp=0)
    return float(a)

def basic_info(G):
    f=open('basic_info.txt','w')
    f.write('网络节点数：')
    f.write(str(G.number_of_nodes()) + '\n')
    f.write('网络边数：')
    f.write(str(G.size()) + '\n')
    f.write('网络边加权和：')
    f.write(str(G.size(weight='weight')) + '\n')
    scc=nx.strongly_connected_components(G)#返回强连通子图的list
    wcc=nx.weakly_connected_components(G)#返回弱连通子图的list
    
  #  f.write('弱连通子图个数：')
  #  f.write(str(len(wcc)) + '\n')
  #  f.write('强连通子图个数：')
  #  f.write(str(len(scc)) + '\n')
  #  largest_scc=scc[0]#返回最大的强连通子图
 #  largest_wcc=wcc[0]
  #  f.write('最大强连通子图节点数：')
  #  f.write(str(len(largest_scc)) + '\n')

  #  f.write('最大弱连通子图节点数：')
  #  f.write(str(largest_wcc.number_of_nodes()) + '\n')
  #  f.write('最大弱连通子图边数：')
   # f.write(str(largest_wcc.size()) + '\n')
  #  nx.write_gexf(largest_wcc,'Gephi_largest_wcc.gexf')

    convert_igraph(G)
    export_graph(G)
 #   print('最大连通子图节点数：'+str(largest_wcc.number_of_nodes())+'\n')
 #   print('最大连通子图边数：'+str(largest_wcc.size())+'\n')

    
##    #f.write('有向图平均路径长度：')
##    #f.write(str(nx.average_shortest_path_length(G)) + '\n')
##    '''
##    f.write('有向图直径：')
##    f.write(str(nx.diameter(G)) + '\n')
##    '''#无法计算有向图直径，因为此时有很多路径长度为无穷大
##    G=G.to_undirected()
##    f.write('平均度：')
##    average_degree=G.size()*2*1.0/G.__len__()
##    f.write(str(average_degree) + '\n')
##    f.write('平均无权聚类系数：')
##    f.write(str(nx.average_clustering(G)) + '\n')#这个命令默认计算无权聚类系数
##    '''
##    f.write('平均无权聚类系数：')
##    f.write(str(nx.average_clustering(G,weight=1)) + '\n')
##    '''
##    #f.write('平均路径长度：')
##    #f.write(str(nx.average_shortest_path_length(G)) + '\n')
##
##    f.write('平均无权相配系数：')
##    f.write(str(nx.degree_assortativity_coefficient(G)) + '\n')
##
##    #f.write('网络直径：')
##    #f.write(str(nx.diameter(G)) + '\n')
##    '''
##    f.write('平均带权相配系数：')
##    f.write(str(nx.degree_assortativity_coefficient(G,weight='weight')) + '\n')
##    '''
##    clustering_co=nx.clustering(G)
##    f=open('clustering_coefficient_average.txt','w')#这个文件输出无权聚类系数的平均值
##    f2=open('clustering_coefficient.txt','w')
##    
##    average_clusterting = {}
##    Times = {}
##    distKeys = range(0,G.__len__()+1)#G.__len__()返回的图中的节点的个数
##    Times = dict([(k,0) for k in distKeys])
##    average_clusterting = dict([(k,0) for k in distKeys])
##    
##    for node in clustering_co:
##        f2.write(str(G.degree(node))+':'+str(clustering_co[node]) + '\n')
##        Times[G.degree(node)] += 1
##        average_clusterting[G.degree(node)] = average_clusterting[G.degree(node)] + clustering_co[node]
##    for each in average_clusterting:
##        if Times[each]!=0:
##            average_clusterting[each] = average_clusterting[each]/Times[each]
##            f.write(str(each)+':'+str(average_clusterting[each]) + '\n')
##
##    #下面生成和G对应的ER Graph
##    f=open('basic_info.txt','a')
##    prob = G.size()*2*1.0 / (1.0*G.number_of_nodes()*(G.number_of_nodes()-1))
##    ER_Graph = nx.random_graphs.erdos_renyi_graph(G.number_of_nodes(),prob)
##    #f.write('ER图平均路径长度：') #由于不是连通图，所以平均最短路径一般无法计算
##    #f.write(str(nx.average_shortest_path_length(ER_Graph)) + '\n')
##    '''
##    虽然ER图不连通，但是可以计算其wcc的聚类系数和路径长度
##    '''
##    print "聚类系数为："+str(nx.average_clustering(ER_Graph)*1.0)
##    f.write('ER图平均无权聚类系数：')
##    f.write(str(nx.average_clustering(ER_Graph)*1.0) + '\n')
##    cc_ER=nx.connected_components(ER_Graph)
##    largest_lcc_ER=ER_Graph.subgraph(cc_ER[0])
##    f.write('ER图的最大连通子图平均无权聚类系数：')
##    f.write(str(nx.average_clustering(largest_lcc_ER)) + '\n')
##    f.write('ER图的最大连通子图平均路径长度：')
##    f.write(str(nx.average_shortest_path_length(largest_lcc_ER)) + '\n')
            

    
    

def cumlutive_degree_distribution(pdf):
    '''
    计算累积度分布
    '''
    pdf_temp=pdf
    scope = range(min(pdf),max(pdf)+1)
    for degree in scope:
        k=degree+1
        while k<=max(pdf):
            pdf[degree]+=pdf_temp[k]
            k+=1
    return pdf

def calling_correlation(G):
    '''
    调用次数和边的相关关系分析
    生成文件的每一条记录中：
    边；边的权重（也就是调用次数）；边的起点的入度；边的终点的入度；边的起点的出度；边的终点的出度。
    '''
    file=open('calling_correlation.txt','w')
    for e in G.edges_iter():
        file.write(str(e)+';')
        '''
        file.write(str(G.get_edge_data(*e)) + ';')
        file.write(str(G.in_degree(e).values()) + ';')
        file.write(str(G.out_degree(e).values()) + '\n')
        '''
        weight=G.get_edge_data(*e)
        for each in weight:
            file.write(str(weight[each]) + ';')
        in_degree=G.in_degree(e)
        for each in  in_degree:
            file.write(str(in_degree[each]) + ';')
        out_degree=G.out_degree(e)
        for each in  out_degree:
            file.write(str(out_degree[each]) + ';')
        file.write('\n')

def calling_correlation_by_node(G):
    '''
    调用次数和节点的相关关系分析
    生成文件的每一条记录中：
    点；点的入度；点的出度；点的所有入边的权重之和（点对应的方法被调用的总次数）；点的所有出边的权重之和（点对应的方法调用别的方法的总次数）
    '''
    file=open('calling_correlation_by_node.txt','w')
    for n in G.nodes_iter():
        file.write(str(n)+';')
        file.write(str(G.in_degree(n))+';')
        file.write(str(G.out_degree(n))+';')
        
        in_edges=G.in_edges(n)
        in_weight=0
        for e in in_edges:
            weight=G.get_edge_data(*e)
            for each in weight:
                in_weight+=weight[each]
        file.write(str(in_weight)+';')
        
        out_edges=G.out_edges(n)
        out_weight=0
        for e in out_edges:
            weight=G.get_edge_data(*e)
            for each in weight:
                out_weight+=weight[each]
        file.write(str(out_weight))
        file.write('\n')

def node_vertex_entropy_compute(G):
    '''
    20130603：节点的熵和点强度记录
    生成文件的每一行记录中：
    点；点的入度；点的入强度（也就是被调用的总次数）；点的入度的entropy；我们自己定义的点的入度的computing_complex；点的整个的local entropy
    '''
    file=open('node_vertex_entropy.txt','w')
    for n in G.nodes_iter():
        file.write(str(n)+';')
	file.write(str(G.in_degree(n))+';')
        
        in_edges=G.in_edges(n)
        in_weight=0
        in_sum=0

        out_edges=G.out_edges(n)
        out_weight=0
        out_sum=0

        vertex_strength=0
        vertex_sum=0
        
        for e in in_edges:
            weight=G.get_edge_data(*e)
            for each in weight:
                in_weight+=weight[each]
                vertex_strength+=weight[each]
				
	file.write(str(in_weight)+';')

        for e in out_edges:
            weight=G.get_edge_data(*e)
            for each in weight:
                out_weight+=weight[each]
                vertex_strength+=weight[each]

        
        for e in in_edges:
            weight=G.get_edge_data(*e)
            for each in weight:
                in_sum+=(weight[each]/in_weight)*(log(weight[each]/in_weight))
                vertex_sum+=(weight[each]/vertex_strength)*(log(weight[each]/vertex_strength))

        for e in out_edges:
            weight=G.get_edge_data(*e)
            for each in weight:
                out_sum+=(weight[each]/out_weight)*(log(weight[each]/out_weight))
                vertex_sum+=(weight[each]/vertex_strength)*(log(weight[each]/vertex_strength))

        if(G.in_degree(n)<=1):
            file.write('0;')
            file.write(str(in_weight)+';')
        else:
            in_entropy=-in_sum/log(G.in_degree(n))
            file.write(str(in_entropy)+';')
            if(G.in_degree(n)>=3):
                adjust=G.in_degree(n)-2
            else:
                adjust=1
            computing_complex=(in_weight*(exp(sqrt(adjust)*in_entropy)))
            file.write(str(computing_complex)+';')
            

        # if(G.out_degree(n)<=1):
            # file.write('0;')
        # else:
            # out_entropy=-out_sum/log(G.out_degree(n))
            # file.write(str(out_entropy)+';')

        if((G.out_degree(n)+G.in_degree(n))==1):
            file.write('0\n')
        else:
            node_entropy=-vertex_sum/log(G.out_degree(n)+G.in_degree(n))
            file.write(str(node_entropy)+'\n')
            
        
        
def degree_analysis(G):
    #print nx.degree_histogram(G)#这里不使用这个函数的原因，是因为其没有区分出度和入度
    #print G.in_degree()
    in_degree = []
    out_degree = []
    list_in=G.in_degree()
    list_out=G.out_degree()
    file_in=open('in_degree.txt','w')
    file_out=open('out_degree.txt','w')
    file_degree=open('degree_correlation.txt','w')

    for n in G.nodes_iter():
        file_degree.write(str(G.in_degree(n)) + ';')
        file_degree.write(str(G.out_degree(n)) + '\n')

    for each_node in list_in:
        file_in.write(each_node+":"+str(list_in[each_node])+"\n")
        in_degree.append(list_in[each_node])

    for each_node in list_out:
        file_out.write(each_node+":"+str(list_out[each_node])+"\n")
        out_degree.append(list_out[each_node])

    
    #print in_degree
    #print out_degree
    
    pdf_in=prob_density_function(in_degree)
    f=open('pdf_in.txt','w')
    for each_degree in pdf_in:
        f.write(str(each_degree)+":"+str(pdf_in[each_degree])+"\n")
    
    pdf_out=prob_density_function(out_degree)
    f=open('pdf_out.txt','w')
    for each_degree in pdf_out:
        f.write(str(each_degree)+":"+str(pdf_out[each_degree])+"\n")
    
    cdf_in=cumlutive_degree_distribution(pdf_in)
    f=open('cdf_in.txt','w')
    for each_degree in cdf_in:
        f.write(str(each_degree)+":"+str(cdf_in[each_degree])+"\n")

    cdf_out=cumlutive_degree_distribution(pdf_out)
    f=open('cdf_out.txt','w')
    for each_degree in cdf_out:
        f.write(str(each_degree)+":"+str(cdf_out[each_degree])+"\n")

    G_Und=G.to_undirected()
    #pdf=nx.degree_histogram(G)#这样得到的是一个list，但是更需要一个显示出每个节点度的list
    #print pdf
    degree_list = []
    degree_dict=G.degree()
    for each_node in degree_dict:
        degree_list.append(degree_dict[each_node])

    pdf=prob_density_function(degree_list)
    f=open('pdf_all.txt','w')
    for each_degree in pdf:
        f.write(str(each_degree)+":"+str(pdf[each_degree])+"\n")

    cdf=cumlutive_degree_distribution(pdf)
    f=open('cdf_all.txt','w')
    for each_degree in cdf:
        f.write(str(each_degree)+":"+str(cdf[each_degree])+"\n")

def I_weighted_triangles_and_degree_iter(G, nodes=None, weight='weight'):
    """ Return an iterator of (node, degree, weighted_triangles).  
    
    Used for weighted clustering.

    """
    if weight is None or G.edges()==[]:
        max_weight=1.0
    else:
        max_weight=float(max(d.get(weight,1.0) 
                             for u,v,d in G.edges(data=True)))
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs= ( (n,G[n]) for n in G.nbunch_iter(nodes) )

    for i,nbrs in nodes_nbrs:
        inbrs=set(nbrs)-set([i])
        weighted_triangles=0.0
        seen=set()
        for j in inbrs:
            wij=G[i][j].get(weight,1.0)
            seen.add(j)
            jnbrs=set(G[j])-seen # this keeps from double counting
            for k in inbrs&jnbrs:
                #wjk=G[j][k].get(weight,1.0)
                wki=G[i][k].get(weight,1.0)
                weighted_triangles+=(wij+wki)/2.0
        yield (i,len(inbrs),weighted_triangles)


def get_dynamic_data(G):
    G_Und=G.to_undirected()
    vs={}#vertex strength
    wc={}#weighted clustering
    #get vertex strength
    for node in G_Und.nodes():
        strength=0.0
        for i in nx.all_neighbors(G_Und,node):
            strength=strength+G_Und[node][i]['weight']
        vs[node]=strength
    #get weighted average Nearest neighbors degree
    test=nx.average_neighbor_degree(G_Und,weight='weight')
    #get unweighted average Nearest neighbors degree
    test2=nx.average_neighbor_degree(G_Und)
    #get unweighted degree correlation function
    test3=nx.average_degree_connectivity(G_Und)
    #get weighted degree correlation function
    test4=nx.average_degree_connectivity(G_Und,weight='weight')
    
    td_iter=I_weighted_triangles_and_degree_iter(G_Und,weight='weight')
    for v,d,t in td_iter:
        if t==0:
            wc[v]=0.0
        else:
            wc[v]=t/(vs[v]*(d-1))
    sac=nx.clustering(G_Und)
    f=open('clustering_coefficient_weighted.txt','w')
    fdegree=open('dgree_correlation_function.txt','w')
    f_cluster_weight_average=open('clustering_coefficient_weighted_average.txt','w')

    average_clusterting = {}
    Times = {}
    distKeys = range(0,G.__len__()+1)#G.__len__()返回的图中的节点的个数
    Times = dict([(k,0) for k in distKeys])
    average_clusterting = dict([(k,0) for k in distKeys])

    total = 0
    total_vertex_strength = 0


    #f.write('点强度  相邻节点度加权平均值      无权度平均值    带权聚集系数   无权聚集系数\n')

    for i in test.keys():
        #f.write(str(i)+','+str(vs[i])+','+str(test[i])+','+str(test2[i])+','+str(wc[i])+','+str(sac[i])+'\n')
        f.write(str(G.degree(i))+':'+str(wc[i])+':'+str(sac[i])+'\n')#wc是带权聚集系数的dict
        Times[G.degree(i)] += 1
        average_clusterting[G.degree(i)] = average_clusterting[G.degree(i)] + wc[i]
        total = total + wc[i]
        total_vertex_strength = total_vertex_strength + vs[i]

    average_weighted_clustering = total / G_Und.__len__()
    average_vertex_strength = total_vertex_strength / G_Und.__len__()
    f_basic_info = open('basic_info.txt','a')#a表示在后面添加，而不是重新写文件
    f_basic_info.write('带权平均聚集系数：'+str(average_weighted_clustering)+'\n')
    f_basic_info.write('平均点强度：'+str(average_vertex_strength)+'\n')

    for each in average_clusterting:
        if Times[each]!=0:
            average_clusterting[each] = average_clusterting[each]/Times[each]
            f_cluster_weight_average.write(str(each)+':'+str(average_clusterting[each]) + '\n')

    

    #fdegree.write('度  加权度相关系数  无权度相关系数\n')
    
    for d in test3.keys():
        fdegree.write(str(d)+','+str(test4[d])+','+str(test3[d])+'\n')
    f.close()
    fdegree.close()
    vertex_strength=[]
    file_vs=open('vertex_strength.txt','w')
    for i in vs:
        file_vs.write(i+":"+str(vs[i])+"\n")
        vertex_strength.append(vs[i])
    
    vertex_strength_int=[int(x) for x in vertex_strength]#使用这条语句可以统一地将float型数组转成int型数组
    #pdf_vertex_strength=prob_density_function(vertex_strength_int)
    #print pdf_vertex_strength
    '''
    cdf_vertex_strength=cumlutive_degree_distribution(pdf_vertex_strength)
    f=open('cdf_vertex_strength.txt','w')
    for each_degree in cdf_vertex_strength:
        f.write(str(each_degree)+":"+str(cdf_vertex_strength[each_degree])+"\n")
    '''
    
    return True


def export_graph(G):
    nx.write_gexf(G,'Gephi.gexf')
    nx.write_pajek(G,'Pajek.net')


def normalized_structure_entropy(G):
    file=open('normalized_structure_entropy.txt','w')
    total_degree=0
    
    for n in G.nodes_iter():
        total_degree=total_degree+G.degree(n)
    
    node_significance=[]
    
    for n in G.nodes_iter():
        if(G.degree(n)==0):
            continue
        node_significance.append(float(G.degree(n))/total_degree)#两个int相除为0，不知道为什么有时候不需要这样进行类型转换？
    structure_entropy=0
    
    #print node_significance
    
    for n in node_significance:
        structure_entropy=structure_entropy-n*log(n)
    #print structure_entropy

    node_number=G.number_of_nodes()
    print "节点数为："+str(node_number)
    normalized_entropy=(structure_entropy-log(4*(node_number-1))/2)/(log(node_number)-log(4*(node_number-1))/2)
    #print normalized_entropy
    file.write(str(normalized_entropy))


def static_analysis(file):
    G = nx.DiGraph()
    findFile = open(sys.argv[1],'r')
    each_lines= findFile.readlines()
    for each_line in each_lines:
        if each_line.__contains__('>'):
            edge=each_line.split('>');
            edge[0]=edge[0][edge[0].index('\"')+1:edge[0].rindex('\"')]
            if(edge[1].__contains__('[')):
                edge[1]=edge[1][edge[1].index('\"'):edge[1].index('[label=')]
                edge[1]=edge[1][edge[1].index('\"')+1:edge[1].rindex('\"')]
                if(G.has_edge(edge[0],edge[1])==False):
                    G.add_edge(edge[0],edge[1])
            else:
                if each_line.count('\"')==2:
                    node=each_line[each_line.index('\"')+1:each_line.rindex('\"')]
                    if(G.has_node(node)==False):
                        G.add_node(node)       
    findFile.close()
    return G

def get_class_name(node):
    method_full_qualifier=node.split(':')
    method_full_name=method_full_qualifier[1].split('.')
    class_name=method_full_name[0]
    return class_name

def community_entropy(communities_quantity):
    quantity_sum=sum(communities_quantity)
    entropy=0
    for each in communities_quantity:
        entropy=entropy-(float(each)/quantity_sum)*log(float(each)/quantity_sum)
    if(len(communities_quantity)!=1):
        entropy=entropy/log(len(communities_quantity))#这里我们可以尝试一下不用归一化的熵值。
    entropy=1-entropy#1减去Entropy以后就是Cohesion Metrics
    #print entropy
    return entropy

def out_put_class_community_results(class_community_result,algorithm):
    file=open('class_fault_cohesion_'+str(algorithm)+'.csv','w')
    file2=open('class_community_result_'+str(algorithm)+'.csv','w')
    for each_class in class_community_result:
        
        file2.write(str(each_class))

        res={}
        community=class_community_result[each_class]
        
        for i in community:
            res[i] = res.get(i, 0) + 1
            file2.write(","+str(i))
        file2.write("\n")

        
        #print each_class
        communities=[k for k in res.keys()]
        communities_quantity=[v for v in res.values()]
        #print max(communities_quantity)
        #print sum(communities_quantity)
        entropy=community_entropy(communities_quantity)
        class_cohesion_old=float(max(communities_quantity))/float(sum(communities_quantity))
        #class_cohesion_old2=float(len(communities))/float(sum(communities_quantity))
        class_cohesion=len(communities)*exp(1-entropy)
        #print class_cohesion
        #if(class_fault_prone.has_key(each_class)):
            #file.write(str(each_class)+','+str(class_fault_prone[each_class])+','+str(class_cohesion)+','+str(entropy)+','+str(len(communities))+"\n")
        file.write(str(each_class)+','+str(sum(communities_quantity))+','+str(class_cohesion_old)+','+str(class_cohesion)+','+str(entropy)+','+str(len(communities))+"\n")
        #目前每一行记录的含义：类名，类的方法数，所在最多的社团的方法数/总的方法数，自己定义的结合社团数+Entropy的内聚度指标，Entropy，所在的社团数

def out_put_community_size(community_size,algorithm):
    f=open('community_size_CDF_'+str(algorithm)+'.txt','w')
    f2=open('community_size_PDF_'+str(algorithm)+'.txt','w')
    pdf_community=prob_density_function(community_size)
    for each_community_size in pdf_community:
        f2.write(str(each_community_size)+":"+str(pdf_community[each_community_size])+"\n")
    
    cdf_community=cumlutive_degree_distribution(pdf_community)

    for each_community_size in cdf_community:
        f.write(str(each_community_size)+":"+str(cdf_community[each_community_size])+"\n")


def correlate_with_class_fault(class_community):
    fault_record=open('class_fault_prone_data.csv','r')
    new_record=open('class_fault_prone_data_new.csv','w')
    class_fault_prone={}
    lines=fault_record.readlines()
    for each_line in lines:
        line_content=each_line.split(',')
        class_name=line_content[0]
        class_name=class_name[class_name.rindex('.')+1:]
        #print class_name
        if(class_community.has_key(class_name)):
            new_record.write(str(class_name)+','+str(line_content[1]))
            fault_status=line_content[1].strip('\n')
            #if(int(fault_status)>0):#看看这样改，能不能用Logistic Regression预测。
                #fault_status=1
            class_fault_prone[class_name]=fault_status

    return class_fault_prone
    

def convert_igraph(G):
    g=igraph.Graph(len(G))
    node_list=G.nodes()
    edge_list=G.edges()
    print('start community analysis')

    #i=1
    #print node_list
    for e in G.edges_iter():
        source=e[0]
        dist=e[1]
        s=node_list.index(source)
        d=node_list.index(dist)
        if(g.are_connected(s,d)==False):#这一条语句是判断g中是否已有这个边，如果不判断，后面的label_propagation算法会报错。
            g.add_edges((s,d))
            #i=i+1
            #print i
    #print g
    print g.summary()
    
##    VertexCluster=g.community_leading_eigenvector()
##    f=open('Community_Result_le.csv','w')
##    f.write("Id,Modularity Class"+"\n")
##    community_index=1
##    for each_community in VertexCluster:
##        print each_community
##        for each_node in each_community:
##            node_name=node_list[each_node]
##            f.write(str(node_name)+","+str(community_index)+"\n")
##        community_index=community_index+1

    print('start community infomap analysis')
    VertexCluster_im=g.community_infomap()
    class_community_result_im={}
    community_size_im=[]
    
    f=open('Community_Result_im.csv','w')
    modularity_score=g.modularity(VertexCluster_im)
    f.write("Id,Modularity Class,Modularity Score:"+str(modularity_score)+",Community Number:"+str(len(VertexCluster_im))+"\n")
    community_index=1
    for each_community in VertexCluster_im:
        #print each_community
        for each_node in each_community:
            node_name=node_list[each_node]
            class_name=get_class_name(node_name)
            if(class_community_result_im.has_key(class_name)):
                class_community_result_im[class_name].append(community_index)
            else:
                class_community_result_im[class_name]=[]
                class_community_result_im[class_name].append(community_index)
            f.write(str(node_name)+","+str(community_index)+"\n")
        community_index=community_index+1
        community_size_im.append(len(each_community))

    #print class_community_result_im
    #print community_size_im
    #class_fault_prone=correlate_with_class_fault(class_community_result_im)
    
    out_put_class_community_results(class_community_result_im,'im')
    #out_put_community_size(community_size_im,'im')


#    print('start community fast greedy analysis')

 
    
    #VertexCluster=g.community_edge_betweenness(None,False,None)
#    VertexDendrogram=g.community_fastgreedy()
#    class_community_result_fg={}
#    community_size_fg=[]
    
##    VertexCluster_fg=VertexDendrogram.as_clustering()
##    modularity_score=g.modularity(VertexCluster_fg)
##    f=open('Community_Result_fg.csv','w')
##    f.write("Id,Modularity Class,Modularity Score:"+str(modularity_score)+",Community Number:"+str(len(VertexCluster_fg))+"\n")
##    community_index=1
##    for each_community in VertexCluster_fg:
        #print each_community
##        for each_node in each_community:
##            node_name=node_list[each_node]
##            class_name=get_class_name(node_name)
##            if(class_community_result_fg.has_key(class_name)):
##                class_community_result_fg[class_name].append(community_index)
##            else:
##                class_community_result_fg[class_name]=[]
##                class_community_result_fg[class_name].append(community_index)
##            f.write(str(node_name)+","+str(community_index)+"\n")
##        community_index=community_index+1
##        community_size_fg.append(len(each_community))

##    out_put_class_community_results(class_community_result_fg,'fg')
    #out_put_community_size(community_size_fg,'fg')


##    print('start community label propagation analysis') 

#    VertexCluster_lp=g.community_label_propagation()
#     class_community_result_lp={}
#     community_size_lp=[]

#     modularity_score=g.modularity(VertexCluster_lp)
#     f=open('Community_Result_lp.csv','w')
#     f.write("Id,Modularity Class,Modularity Score:"+str(modularity_score)+",Community Number:"+str(len(VertexCluster_lp))+"\n")
#     community_index=1
#     for each_community in VertexCluster_lp:
        #print each_community
#         for each_node in each_community:
#             node_name=node_list[each_node]
#             class_name=get_class_name(node_name)
#             if(class_community_result_lp.has_key(class_name)):
#                 class_community_result_lp[class_name].append(community_index)
#             else:
 #                class_community_result_lp[class_name]=[]
#                 class_community_result_lp[class_name].append(community_index)
            
#             f.write(str(node_name)+","+str(community_index)+"\n")
 #        community_index=community_index+1
#         community_size_lp.append(len(each_community))

 #    out_put_class_community_results(class_community_result_lp,'lp')
    #out_put_community_size(community_size_lp,'lp')

    #print igraph.compare_communities(VertexCluster_im,VertexCluster_fg,'nmi',True)

 #    print('start community multilevel analysis')
  #   VertexCluster_ml=g.community_multilevel()
 #    class_community_result_ml={}
#     community_size_ml=[]

 #    modularity_score=g.modularity(VertexCluster_ml)
 #    f=open('Community_Result_ml.csv','w')
 #    f.write("Id,Modularity Class,Modularity Score:"+str(modularity_score)+",Community Number:"+str(len(VertexCluster_ml))+"\n")
 #    community_index=1
#     for each_community in VertexCluster_ml:
        #print each_community
  #       for each_node in each_community:
 #            node_name=node_list[each_node]
 #            class_name=get_class_name(node_name)
 #            if(class_community_result_ml.has_key(class_name)):
 #                class_community_result_ml[class_name].append(community_index)
 #            else:
 #                class_community_result_ml[class_name]=[]
 #                class_community_result_ml[class_name].append(community_index)
 #            
 #            f.write(str(node_name)+","+str(community_index)+"\n")
#         community_index=community_index+1
#         community_size_ml.append(len(each_community))

#     out_put_class_community_results(class_community_result_ml,'ml')
    #out_put_community_size(community_size_ml,'ml')

 #    f=open('normalized_mutual_information_results.csv','w')
 #    f.write(str('fg-im,'+str(igraph.compare_communities(VertexCluster_im,VertexCluster_fg,'nmi')))+"\n")
 #    f.write(str('fg-lp,'+str(igraph.compare_communities(VertexCluster_lp,VertexCluster_fg,'nmi')))+"\n")
 #    f.write(str('fg-ml,'+str(igraph.compare_communities(VertexCluster_ml,VertexCluster_fg,'nmi')))+"\n")
#     f.write(str('im-lp,'+str(igraph.compare_communities(VertexCluster_im,VertexCluster_lp,'nmi')))+"\n")
 #    f.write(str('im-ml,'+str(igraph.compare_communities(VertexCluster_im,VertexCluster_ml,'nmi')))+"\n")
 #    f.write(str('lp-ml,'+str(igraph.compare_communities(VertexCluster_lp,VertexCluster_ml,'nmi')))+"\n")

##    VertexCluster=g.community_edge_betweenness()
##    modularity_score=g.modularity(VertexCluster)
##    f=open('Community_Result_eb.csv','w')
##    f.write("Id,Modularity Class,Modularity Score:"+str(modularity_score)+",Community Number:"+str(len(VertexCluster))+"\n")
##    community_index=1
##    for each_community in VertexCluster:
##        #print each_community
##        for each_node in each_community:
##            node_name=node_list[each_node]
##            f.write(str(node_name)+","+str(community_index)+"\n")
##        community_index=community_index+1 
    


file = sys.argv[1]
G = static_analysis(file)
basic_info(G)
##degree_analysis(G)
##calling_correlation(G)
##calling_correlation_by_node(G)
###get_dynamic_data(G)
##node_vertex_entropy_compute(G)
##
###export_graph(G)
##normalized_structure_entropy(G)
#convert_igraph(G)
