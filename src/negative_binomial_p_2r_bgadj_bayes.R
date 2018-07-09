library(LSD)
### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_track_file = args[1]
signal_folder = args[2]
input_track_file = args[3]
input_folder = args[4]
output_name = args[5]

mean_vec = c()
var_vec = c()
size_vec = c()
prob_vec = c()

get_true_NB_prob_size = function(mean_non0, mean_x2_non0){
	###### identify p0
	best_p0 = 0
	best_prob_dif = 1
	k=0
	for (i in seq(0,0.99,0.005)){
		k = k+1
		ProbT = mean_non0 / (mean_x2_non0 - mean_non0^2 * (1-i))
		SizeT = mean_non0 * (1-i) * ProbT / (1-ProbT)
		nb_v_T = rnbinom(1e+4, SizeT, ProbT)
		p0_new = sum(nb_v_T==0) / length(nb_v_T)
		p0_dif = abs(i-p0_new)
		if ((k%%100)==0){
			print(paste('iteration:', toString(k)))
		}
		if (abs(i-p0_new) < best_prob_dif){
			print(paste('iteration:', toString(k)))
			print('change best_p0')
			best_prob_dif = abs(i-p0_new)
			best_p0 = i
		}
	}
	print('true p0:')
	print(sum(nb_v==0)/length(nb_v))
	print('estimated p0: ')
	print(best_p0)
	ProbT = mean_non0 / (mean_x2_non0 - mean_non0^2 * (1-best_p0))
	SizeT = mean_non0 * (1-best_p0) * ProbT / (1-ProbT)
	return(c(ProbT, SizeT, best_p0))
}

### read data
sig = read.table(paste(signal_folder, signal_track_file, sep=''), header = F)
input = read.table(paste(input_folder, input_track_file, sep=''), header = F)
thesh = 0

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
### get sig bg regions no bgs
sig_bg = sig[,1]
sig_bg_non0 = sig_bg[sig_bg>thesh]
sig_bg_mean = mean(sig_bg_non0)
sig_bg_mean_sig2 = mean(sig_bg_non0^2)
sig_bg_var = var(sig_bg_non0)

probT_sizeT_p0 = get_true_NB_prob_size(sig_bg_mean, sig_bg_mean_sig2)

print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_bg_var/sig_bg_mean, digits=3)) ))
print(sig_bg_mean)
print(sig_bg_var)
print(length(sig_bg_non0))

### get negative binomial parameters from signal track bg regions
sig_bg_prob = probT_sizeT[1]
if (sig_bg_prob<0.1){
	sig_bg_prob = 0.1
}

if (sig_bg_prob>=0.9){
	sig_bg_prob = 0.9
}

p0 = probT_sizeT[3]
sig_bg_size = sig_bg_mean^2 * (1-p0) / (sig_bg_mean_sig2 - sig_bg_mean^2 * (1-p0) - sig_bg_mean)
### get input bg regions
input_bg = input[,1]
input_bg_non0 = input_bg[input_bg>thesh]
input_bg_mean = mean(input_bg_non0)
inpy_bg_var = var(input_bg_non0)
print(paste('check input track overdispersion in background regions, var/mean=', toString(round(inpy_bg_var/input_bg_mean, digits=3)) ))
print(sig_bg_prob)
print(sig_bg_size)
print(length(input_bg_non0))

print(head(input_bg))
print(summary(input_bg))
print(input_bg_mean)
print(inpy_bg_var)

### get negative binomial p-value
sig_input = cbind(sig, input)
nb_pval = apply(sig_input, MARGIN=1, function(x) pnbinom(x[1], sig_bg_size, sig_bg_prob, lower.tail=FALSE) )
### get -log10(p-value)
print('get -log10(p-value)')
print(min(nb_pval[nb_pval!=0]))
print(length(nb_pval))
print(length(nb_pval[nb_pval<=1e-100]))
nb_pval_min = min(nb_pval[nb_pval!=0])
nb_pval[nb_pval<=1e-100] = 1e-100
print(summary(nb_pval))
############### second round
### get sig bg regions
sig_bg = sig[nb_pval>=0.001,]
sig_bg_non0 = sig_bg[sig_bg>thesh]
sig_bg_mean = mean(sig_bg_non0)
sig_bg_mean_sig2 = mean(sig_bg_non0^2)
sig_bg_var = var(sig_bg_non0)

probT_sizeT_p0 = get_true_NB_prob_size(sig_bg_mean, sig_bg_mean_sig2)

print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_bg_var/sig_bg_mean, digits=3)) ))
print(sig_bg_mean)
print(sig_bg_var)
print(length(sig_bg_non0))

### get negative binomial parameters from signal track bg regions
sig_bg_prob = probT_sizeT[1]
if (sig_bg_prob<0.1){
	sig_bg_prob = 0.1
}

if (sig_bg_prob>=0.9){
	sig_bg_prob = 0.9
}

p0 = probT_sizeT[3]
sig_bg_size = sig_bg_mean^2 * (1-p0) / (sig_bg_mean_sig2 - sig_bg_mean^2 * (1-p0) - sig_bg_mean)

mean_vec[1] = sig_bg_mean
var_vec[1] = sig_bg_var
size_vec[1] = sig_bg_size
prob_vec[1] = sig_bg_prob

### get input bg regions
input_bg = input[nb_pval>=0.001,]
input_bg_non0 = input_bg[input_bg>thesh]
input_bg_mean = mean(input_bg_non0)
inpy_bg_var = var(input_bg_non0)
print(paste('check input track overdispersion in background regions, var/mean=', toString(round(inpy_bg_var/input_bg_mean, digits=3)) ))
print(sig_bg_prob)
print(sig_bg_size)
print(length(input_bg_non0))

print(head(input_bg))
print(summary(input_bg))
print(input_bg_mean)
print(inpy_bg_var)

### get negative binomial p-value
sig_input = cbind(sig, input)
nb_pval = apply(sig_input, MARGIN=1, function(x) pnbinom(x[1], sig_bg_size * (x[2]+1)/(input_bg_mean+1), sig_bg_prob, lower.tail=FALSE) )
### get -log10(p-value)
nb_pval[nb_pval<=1e-100] = 1e-100
neglog10_nb_pval = -log10(nb_pval)

### write output
write.table(neglog10_nb_pval, paste(output_name, '.nbp_2r_bgadj.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################

info_matrix = cbind(mean_vec, var_vec, size_vec, prob_vec)
write.table(info_matrix, paste(output_name, '.mvsp.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

