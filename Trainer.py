import numpy as np
from The_Payne import training
from The_Payne import utils
import matplotlib.pyplot as plt
# load the default training set. Note that, due to the GitHub size limit,
# this training set is a small subset of what I used to train the default network 
# training_labels1, training_spectra1, validation_labels1, validation_spectra1 = utils.load_training_data()
training_data=np.load('Training_data_no_degredation_all_clean_Fe_Alpha.npy',allow_pickle=True)
#validation=np.load('Validation_data_no_degredation_all_clean.npy',allow_pickle=True)
print(len(training_data))
#training_data=training_data[:int(len(training_data)/2)]
print(len(training_data))
# label array unit = [n_spectra, n_labels]
# spectra_array unit = [n_spectra, n_pixels]

# The validation set is used to independently evaluate how well the neural net
# is emulating the spectra. If the network overfits the spectral variation, while 
# the loss will continue to improve for the training set, the validation set
# should exhibit a worsen loss.

# the codes outputs a numpy array ""NN_normalized_spectra.npz" 
# which stores the trained network parameters
# and can be used to substitute the default one in the directory neural_nets/
# it will also output a numpy array "training_loss.npz"
# which stores the progression of the training and validation losses
training_labels=np.vstack(training_data[:,1])
training_labels=np.hstack((training_labels[:,0:5],training_labels[:,6:]))
training_spectra=np.vstack(np.vstack(training_data[:,0]))
# training_spectra=np.vstack
# validation_labels=np.vstack(validation[:,1])
# validation_labels=np.hstack((validation_labels[:,0:5],validation_labels[:,6:]))
# validation_spectra=np.vstack(np.vstack(validation[:,0]))

training_labels=training_labels.astype('float64')
# validation_labels=validation_labels.astype('float64')

# training_spectra=training_spectra.astype('float64')
# validation_spectra=validation_spectra.astype('float64')

training_spectra=[x[3].astype(float) for x in training_spectra]
training_spectra=np.vstack(training_spectra)

validation_spectra=training_spectra[-len(training_spectra)//4:]
validation_labels=training_labels[-len(training_labels)//4:]

training_spectra=training_spectra[:-len(training_spectra)//4]
training_labels=training_labels[:-len(training_labels)//4]


training.neural_net(training_labels, training_spectra,\
                    validation_labels, validation_spectra,\
                    num_neurons=300, learning_rate=1e-4,\
                    num_steps=1.0e7, batch_size=512,num_pixel=22080,num_features=64*5)

# a larger batch_size (e.g. 512) when possible is desirable
# here we choose batch_size=128 above because the sample training set is limited in size

tmp = np.load("training_loss.npz") # the output array also stores the training and validation loss
training_loss = tmp["training_loss"]
validation_loss = tmp["validation_loss"]

plt.figure(figsize=(14, 4))
plt.plot(np.arange(training_loss.size)*100, training_loss, 'k', lw=0.5, label = 'Training set')
plt.plot(np.arange(training_loss.size)*100, validation_loss, 'r', lw=0.5, label = 'Validation set')
plt.legend(loc = 'best', frameon = False, fontsize= 18)
plt.yscale('log')
# plt.ylim([5,100])
plt.xlabel("Step", size=20)
plt.ylabel("Loss", size=20)
