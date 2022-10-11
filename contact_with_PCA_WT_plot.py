import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import matplotlib
matplotlib.use('TkAgg')
import seaborn as sns 
sns.set_context("talk")
# plt.style.use("dark_background")

def contact_count(contactList,bins,dataName):
    cv_edges = np.histogram_bin_edges(contactList['pca1'],bins=bins,weights=None)
    # print(cv_edges)
    data_list = []
    for i in range(cv_edges.size-1):
        cv_list_bin = contactList.query('%s<=pca1<%s'%(cv_edges[i],cv_edges[i+1]))
        group1_count = cv_list_bin.query('group==1')
        group2_count = cv_list_bin.query('group==2')
        group3_count = cv_list_bin.query('group==3')
        # print(cv_list_bin)

        normalized_factor = len(cv_list_bin.drop_duplicates('frame'))
        # print("%s<=pca1<%s"%(cv_edges[i],cv_edges[i+1]))
        # # print(cv_list_bin)
        # print(cv_list_bin.drop_duplicates('frame'))
        # print("N.F. =  %s"%(normalized_factor))

        if normalized_factor ==0:
            data_list.append((cv_edges[i],0,0,0,dataName))
        else:
            data_list.append((cv_edges[i],len(group1_count)/normalized_factor,
                              len(group2_count)/normalized_factor,len(group3_count)/normalized_factor,dataName))
    data_list = np.stack(data_list)
    return data_list

path=''

data_name = 'WT_mwMETAD_contact_count'
system_list = ["WT"]

contacts1 = pd.read_csv(path+data_name+'.dat',delim_whitespace=True)

data = contact_count(contacts1,25,system_list[0])

data = pd.DataFrame(data=data,columns=['pca1','g1','g2','g3','system'])
data.to_csv(path+data_name+'_plot.dat',sep='\t')

######
plot_data = pd.read_csv(path+data_name+'_plot.dat',delim_whitespace='True')

# print(plot_data)
color_theme = sns.color_palette("Dark2")
color_list = ['#F26622','#055053','#EFA320']
linewidth = 4
y_axis_list = ['g1','g2','g3']
marker_list = ['o','^','v','s']

fig,ax_list = plt.subplots()

for a in range(3):
    g = sns.lineplot(x='pca1',y='%s'%(y_axis_list[a]),data=plot_data[['pca1','%s'%(y_axis_list[a])]],
             linewidth=linewidth,color=color_list[a],label='%s'%(y_axis_list[a]),
             marker='X',markersize=15)

    # g.set_xlim([-41,13])
    # g.set_ylim([0,6])
    g.invert_xaxis()
    g.set_ylabel("Average Number of Contacts")
    g.set_xlabel("pca1")
    g.set_title(data_name)

# plt.legend()
plt.rcParams['figure.constrained_layout.use']=True
plt.rcParams['ps.useafm']=True
plt.rcParams['pdf.fonttype']=42
fig.set_size_inches(7.5,6)
plt.locator_params(axis='y',nbins=4)

# plt.savefig(path+str(data_name)+'_plot.png',dpi=300,transparent=False)
# plt.savefig(path+str(data_name)+'_plot.eps',dpi=300,transparent=False)

plt.show()
