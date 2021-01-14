datasets = csvread('./tmp_data.csv');
[beta,support,weight,cons]=aftmix(datasets);
csvwrite('./tmp_result.csv', beta')
csvwrite('./tmp_result_mixing.csv', [support;weight])
csvwrite('./tmp_result_cons.csv', cons)


