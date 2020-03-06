using MAT, DataFrames

function load_data()

	vars = matread("../Data/SGU/sgu_tot_data.mat");

	mat = vars["datamatrix"]
	countrynames = vars["country_name"]

	data = zeros(size(mat,1)*size(mat,3), size(mat,2)+1)
	for jj in 1:size(mat,2)

		var_mat = [mat[jt, jj, jc] for jt in 1:size(mat,1), jc in 1:size(mat,3)]
		data[:, jj] = vec(var_mat)
	end

	varnames = ["year", "tot", "tb", "y", "c", "i", "rer", "BAAFFspread", "country"]

	country_mat = [jc for jt in 1:size(mat,1), jc in 1:size(mat,3)]

	data[:, end] = vec(country_mat)

	df = DataFrame(year = data[:,1], tot = data[:,2], tb = data[:,3], y = data[:,4], c = data[:,5], i = data[:,6], rer = data[:,7], BAAFFspread = data[:,8], country = countrynames[Int.(data[:,9])])

	return df, countrynames
end