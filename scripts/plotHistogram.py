"""
Small script made to plot histograms of ORF lengths predicted by glimmer 
"""

# import matplotlib for making the histogram. 
import matplotlib.pyplot as plt

''' 
TODO:  Change the names of the input and output files here, so that the script uses your glimmer predictions as input,
and outputs a histogram named according to the input genome predictions.
'''
glimmer_pred = "17.glimmer.predict"
out_file = "ORFlength_17_hist.png"

# start an empty list that we can store genelengths in
orf_lengths = list()

# read file and loop over each line
with open(glimmer_pred, "r") as fil:
    # skip the first line to avoid the >name line
    next(fil)
    for line in fil:
        # split the line on whitespace. Note that glimmer returns a multispace seperated file not "\t" or ","
        # Note that python str.split() will split on and remove all whitespace from a string 
        Sline = line.split()

        # calculate the length of the predicted orf
        orf_len = float(Sline[2])-float(Sline[1])
        # Make the length absolute. Since some orfs are on the reverse strand, we do not want negative values.
        orf_len = abs(orf_len)
        # append the orf length to the list of lengths 
        orf_lengths.append(orf_len)

'''
TODO: use matplotlib to plot the orf_lengths in a histogram
'''
plt.hist(orf_lengths,bins=[0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,8000,10000,12000])
#plt.hist(orf_lengths)
plt.xlabel('Predicted gene length (bp)')
plt.ylabel('Frequency')
plt.title("Distribution of predicted gene length (17.fa)")

'''
TODO: Save the histogram to a png file named as the out_file variable
'''
plt.savefig(out_file)
plt.show()
#orf_lengths.sort()
#print(orf_lengths)

