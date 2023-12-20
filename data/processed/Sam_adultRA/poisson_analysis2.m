function [q_maxs_acoustic_groupings2, H_maxs_acoustic_groupings2 ] = poisson_analysis2(case_id)

load(case_id, 'spiketrains', 'number_of_trials', 'ordering', 'groupings','acoustic_data');
number_of_categories = length(groupings)/2;

c = spiketrains{1};

n_random = 2;
q_maxs_acoustic_groupings2 = zeros(n_random,4);
H_maxs_acoustic_groupings2 = zeros(n_random,4);

for count = 1:n_random
    
    r = get_random_spiketrains(c);
    
    distance_matrixes_r = get_distance_matrixes(r);
    distance_matrixes_r_norm = normalize_distance_matrixes(distance_matrixes_r,r);
    
    [N_allzs_r_3dacoustics, H_allzs_r_norm_3dacoustics] = clustering_allzs(distance_matrixes_r_norm, ordering, groupings);
    [N_allzs_r_pitch, H_allzs_r_norm_pitch] = clustering_allzs(distance_matrixes_r_norm, acoustic_data{1}{2}, even_groupings(number_of_trials,number_of_categories));
    [N_allzs_r_amplitude, H_allzs_r_norm_amplitude] = clustering_allzs(distance_matrixes_r_norm, acoustic_data{2}{2}, even_groupings(number_of_trials,number_of_categories));
    [N_allzs_r_spectral_entropy, H_allzs_r_norm_spectral_entropy] = clustering_allzs(distance_matrixes_r_norm, acoustic_data{3}{2}, even_groupings(number_of_trials,number_of_categories));
    
    
    %H_allz_randoms for significance test
    n_random_toadd = 500;

    H_allzs_randoms_r_norm_3dacoustics = sig_test_random(n_random_toadd,distance_matrixes_r_norm,groupings);
    H_allzs_randoms_r_norm_pitch = sig_test_random(n_random_toadd,distance_matrixes_r_norm,even_groupings(number_of_trials,number_of_categories));


    
    
    for_3dacoustics = 1;
        H_allzs_r_norm = H_allzs_r_norm_3dacoustics;
        H_allzs_randoms_r_norm = H_allzs_randoms_r_norm_3dacoustics;
    
        H_random_by_zs = separate_H_randoms(H_allzs_randoms_r_norm);

        H_allzs = H_allzs_r_norm;
        H_allzs_randoms = H_allzs_randoms_r_norm;

        n_random = size(H_allzs_randoms,1);
        H_random = zeros(n_random,11);

        q_values = [0,.05,.1,.2,.3,.5,1,2,5,10,20];
        q_indexes = [1 2 3 4 5 6 7 8 9 10 11];
        
        [H_toplot,max_z_index] = max(H_allzs,[],1);
        
        for i = 1:n_random
            H_allzs_random = H_allzs_randoms{i};
            H_random(i,:)=[H_allzs_random(max_z_index(1),1),H_allzs_random(max_z_index(2),2),H_allzs_random(max_z_index(3),3),H_allzs_random(max_z_index(4),4),H_allzs_random(max_z_index(5),5),H_allzs_random(max_z_index(6),6),H_allzs_random(max_z_index(7),7),H_allzs_random(max_z_index(8),8),H_allzs_random(max_z_index(9),9),H_allzs_random(max_z_index(10),10),H_allzs_random(max_z_index(11),11)];    
        end 
        
        means = mean(H_random,1);
        H_random_ordered = sort(H_random,1);
        H_toplot_95 = H_random_ordered(round(.95*n_random),:); 


        subtractedmean = H_toplot-means;
        subtractedmean(subtractedmean < 0) =0;
        subtractedmean_95 = H_toplot_95-means;

        H_mean_subtracted_3dacoustics = subtractedmean
        significance_of_H_3dacoustics = subtractedmean > subtractedmean_95
        
        
    for_3dacoustics= 0; 
        H_allzs_r_norm = H_allzs_r_norm_pitch;
        H_allzs_randoms_r_norm = H_allzs_randoms_r_norm_pitch;
        
        H_random_by_zs = separate_H_randoms(H_allzs_randoms_r_norm);

        H_allzs = H_allzs_r_norm;
        H_allzs_randoms = H_allzs_randoms_r_norm;

        n_random = size(H_allzs_randoms,1);
        H_random = zeros(n_random,11);

        q_values = [0,.05,.1,.2,.3,.5,1,2,5,10,20];
        q_indexes = [1 2 3 4 5 6 7 8 9 10 11];
        
        [H_toplot,max_z_index] = max(H_allzs,[],1);
        
        for i = 1:n_random
            H_allzs_random = H_allzs_randoms{i};
            H_random(i,:)=[H_allzs_random(max_z_index(1),1),H_allzs_random(max_z_index(2),2),H_allzs_random(max_z_index(3),3),H_allzs_random(max_z_index(4),4),H_allzs_random(max_z_index(5),5),H_allzs_random(max_z_index(6),6),H_allzs_random(max_z_index(7),7),H_allzs_random(max_z_index(8),8),H_allzs_random(max_z_index(9),9),H_allzs_random(max_z_index(10),10),H_allzs_random(max_z_index(11),11)];    
        end 
        
        means = mean(H_random,1);
        H_random_ordered = sort(H_random,1);
        H_toplot_95 = H_random_ordered(round(.95*n_random),:); 


        subtractedmean = H_toplot-means;
        subtractedmean(subtractedmean < 0) =0;
        subtractedmean_95 = H_toplot_95-means;

        H_mean_subtracted_pitch = subtractedmean
        significance_of_H_pitch = subtractedmean > subtractedmean_95
        
    for_3dacoustics = 2; 
        H_allzs_r_norm = H_allzs_r_norm_amplitude;
        H_allzs_randoms_r_norm = H_allzs_randoms_r_norm_pitch; 
        
        H_random_by_zs = separate_H_randoms(H_allzs_randoms_r_norm);

        H_allzs = H_allzs_r_norm;
        H_allzs_randoms = H_allzs_randoms_r_norm;

        n_random = size(H_allzs_randoms,1);
        H_random = zeros(n_random,11);

        q_values = [0,.05,.1,.2,.3,.5,1,2,5,10,20];
        q_indexes = [1 2 3 4 5 6 7 8 9 10 11];
        
        [H_toplot,max_z_index] = max(H_allzs,[],1);
        
        for i = 1:n_random
            H_allzs_random = H_allzs_randoms{i};
            H_random(i,:)=[H_allzs_random(max_z_index(1),1),H_allzs_random(max_z_index(2),2),H_allzs_random(max_z_index(3),3),H_allzs_random(max_z_index(4),4),H_allzs_random(max_z_index(5),5),H_allzs_random(max_z_index(6),6),H_allzs_random(max_z_index(7),7),H_allzs_random(max_z_index(8),8),H_allzs_random(max_z_index(9),9),H_allzs_random(max_z_index(10),10),H_allzs_random(max_z_index(11),11)];    
        end 
        
        means = mean(H_random,1);
        H_random_ordered = sort(H_random,1);
        H_toplot_95 = H_random_ordered(round(.95*n_random),:); 


        subtractedmean = H_toplot-means;
        subtractedmean(subtractedmean < 0) =0;
        subtractedmean_95 = H_toplot_95-means;

        H_mean_subtracted_amplitude = subtractedmean;
        significance_of_H_amplitude = subtractedmean > subtractedmean_95;
        
        
    for_3dacoustics = 3;
        H_allzs_r_norm = H_allzs_r_norm_spectral_entropy;
        H_allzs_randoms_r_norm = H_allzs_randoms_r_norm_pitch;
        
        H_random_by_zs = separate_H_randoms(H_allzs_randoms_r_norm);

        H_allzs = H_allzs_r_norm;
        H_allzs_randoms = H_allzs_randoms_r_norm;

        n_random = size(H_allzs_randoms,1);
        H_random = zeros(n_random,11);

        q_values = [0,.05,.1,.2,.3,.5,1,2,5,10,20];
        q_indexes = [1 2 3 4 5 6 7 8 9 10 11];
        
        [H_toplot,max_z_index] = max(H_allzs,[],1);
        
        for i = 1:n_random
            H_allzs_random = H_allzs_randoms{i};
            H_random(i,:)=[H_allzs_random(max_z_index(1),1),H_allzs_random(max_z_index(2),2),H_allzs_random(max_z_index(3),3),H_allzs_random(max_z_index(4),4),H_allzs_random(max_z_index(5),5),H_allzs_random(max_z_index(6),6),H_allzs_random(max_z_index(7),7),H_allzs_random(max_z_index(8),8),H_allzs_random(max_z_index(9),9),H_allzs_random(max_z_index(10),10),H_allzs_random(max_z_index(11),11)];    
        end 
        
        means = mean(H_random,1);
        H_random_ordered = sort(H_random,1);
        H_toplot_95 = H_random_ordered(round(.95*n_random),:); 


        subtractedmean = H_toplot-means;
        subtractedmean(subtractedmean < 0) =0;
        subtractedmean_95 = H_toplot_95-means;

        H_mean_subtracted_spectral_entropy = subtractedmean
        significance_of_H_spectral_entropy = subtractedmean > subtractedmean_95

    
        if(length(significance_of_H_3dacoustics)>1 && sum(significance_of_H_3dacoustics(1:11)>=1))
            [H_max_3dacoustics, index] = max(H_mean_subtracted_3dacoustics(significance_of_H_3dacoustics))
            q_indexes2 = q_indexes(significance_of_H_3dacoustics);
            q_max_3dacoustics = q_values(q_indexes2(index));
        
        else
            q_max_3dacoustics = -1;
            H_max_3dacoustics = -1;
        end
        
         if(length(significance_of_H_pitch)>1 && sum(significance_of_H_pitch(1:11)>=1))
            [H_max_pitch, index] = max(H_mean_subtracted_pitch(significance_of_H_pitch));
            q_indexes2 = q_indexes(significance_of_H_pitch);
            q_max_pitch = q_values(q_indexes2(index));
           
         else
            q_max_pitch = -1;
            H_max_pitch = -1;
         end
         
         if(length(significance_of_H_amplitude)>1 && sum(significance_of_H_amplitude(1:11)>=1))
            [H_max_amplitude, index] = max(H_mean_subtracted_amplitude(significance_of_H_amplitude))
            q_indexes2 = q_indexes(significance_of_H_amplitude);
            q_max_amplitude = q_values(q_indexes2(index));
        
        else
            q_max_amplitude = -1;
            H_max_amplitude = -1;
        end
        
         if(length(significance_of_H_spectral_entropy)>1 && sum(significance_of_H_spectral_entropy(1:11)>=1))
            [H_max_spectral_entropy, index] = max(H_mean_subtracted_spectral_entropy(significance_of_H_spectral_entropy));
            q_indexes2 = q_indexes(significance_of_H_spectral_entropy);
            q_max_spectral_entropy = q_values(q_indexes2(index))
           
         else
            q_max_spectral_entropy = -1;
            H_max_spectral_entropy = -1;
         end

         q_maxs_acoustic_groupings2(count,:) = [q_max_3dacoustics q_max_pitch q_max_amplitude q_max_spectral_entropy]
         H_maxs_acoustic_groupings2(count,:) = [H_max_3dacoustics H_max_pitch H_max_amplitude H_max_spectral_entropy]
         
end


end



function [r] = get_random_spiketrains(c)

r = cell(size(c));


for i = 1:length(c)
    num_spikes = length(c{i});

    if(num_spikes==0)
        random_spiketrain = [];
    else
        windows = 2;
        while(max(windows)>1)
            random_spiketrain = 40.*rand(1,num_spikes);
            random_spiketrain = sort(random_spiketrain, 'ascend');
            windows = get_windows({random_spiketrain},1);
        end
    end

    r{i} = random_spiketrain;

end
end

function [distance_matrixes] = get_distance_matrixes(e)



    
    n = length(e);


    distance_matrix_q0=zeros(n,n);
    distance_matrix_q1=zeros(n,n);

    distance_matrix_qp1=zeros(n,n);
    distance_matrix_qp3=zeros(n,n);

    distance_matrix_q2=zeros(n,n);
    distance_matrix_q5=zeros(n,n);
    distance_matrix_q10=zeros(n,n);

    distance_matrix_qp2=zeros(n,n);
    distance_matrix_qp5=zeros(n,n);
    distance_matrix_qp05=zeros(n,n);
    distance_matrix_q20=zeros(n,n);

    
        for i = 1:n
            i;
            for j = i+1:n
  
            distance_matrix_q0(i,j) = abs(length(e{i})-length(e{j}));
            distance_matrix_q1(i,j) = VPcost(e{i},e{j},1);

            distance_matrix_qp1(i,j) = VPcost(e{i},e{j},.1);
            distance_matrix_qp3(i,j) = VPcost(e{i},e{j},.3);

            distance_matrix_q2(i,j)=VPcost(e{i},e{j},2);
            distance_matrix_q5(i,j)=VPcost(e{i},e{j},5);
            distance_matrix_q10(i,j)=VPcost(e{i},e{j},10);

            distance_matrix_qp2(i,j) = VPcost(e{i},e{j},.2);
            distance_matrix_qp5(i,j) = VPcost(e{i},e{j},.5);
            distance_matrix_qp05(i,j) = VPcost(e{i},e{j},.05);
            distance_matrix_q20(i,j)=VPcost(e{i},e{j},20);



            distance_matrix_q0(j,i) = distance_matrix_q0(i,j);
            distance_matrix_q1(j,i) = distance_matrix_q1(i,j);
            distance_matrix_qp1(j,i) = distance_matrix_qp1(i,j);
            distance_matrix_qp3(j,i) = distance_matrix_qp3(i,j);
            distance_matrix_qp2(j,i) = distance_matrix_qp2(i,j);
            distance_matrix_qp5(j,i) = distance_matrix_qp5(i,j);
            distance_matrix_q2(j,i) = distance_matrix_q2(i,j);
            distance_matrix_q5(j,i) = distance_matrix_q5(i,j);
            distance_matrix_q10(j,i) = distance_matrix_q10(i,j);
            distance_matrix_q20(j,i) = distance_matrix_q20(i,j);
            distance_matrix_qp05(j,i) = distance_matrix_qp05(i,j);

            end
        end
     

        
        distance_matrixes = {distance_matrix_q0;
                                distance_matrix_qp05;
                                distance_matrix_qp1;
                                distance_matrix_qp2;
                                distance_matrix_qp3;
                                distance_matrix_qp5;
                                distance_matrix_q1;
                                distance_matrix_q2;
                                distance_matrix_q5;
                                distance_matrix_q10;
                                distance_matrix_q20;}
                            

                            
      


end

function [distance_matrixes_norm] = normalize_distance_matrixes(distance_matrixes,d)




q_strings = {'0','p05','p1','p2','p3','p5','1','2','5','10','20'};


[n,n] = size(distance_matrixes{1});



for q_index = 1:length(q_strings)

    distance_matrix = distance_matrixes{q_index};

    distance_matrix_norm = zeros(n,n); 
    
    
    for i = 1:n
        for j = 1:n
            if (length(d{i})+length(d{j})~=0)
                distance_matrix_norm(i,j) = distance_matrix(i,j)/(length(d{i}) + length(d{j}));
            else
                distance_matrix_norm(i,j) = 0;
            end
        end
    end
     eval(char(['distance_matrix_q' q_strings{q_index} ' = distance_matrix_norm;']));       

end

      distance_matrixes_norm = {distance_matrix_q0;
                                distance_matrix_qp05;
                                distance_matrix_qp1;
                                distance_matrix_qp2;
                                distance_matrix_qp3;
                                distance_matrix_qp5;
                                distance_matrix_q1;
                                distance_matrix_q2;
                                distance_matrix_q5;
                                distance_matrix_q10;
                                distance_matrix_q20;}
      
             
  

end


function [N_allzs,H_allzs] = clustering_allzs(distance_matrixes,ordering,groupings)
    % [N_allzs, H_allzs] is output 

    z_values = [-8,-4,-2,-1,1,2,4,8];
    
    N_allzs = cell(8,1);
    H_allzs = zeros(8,11);
    
    

    for z_index = 1:8
        [N,H]=clustering_withdistancematrixes(distance_matrixes,ordering, z_values(z_index),0, groupings);
        N_allzs{z_index,:} = N;
        H_allzs(z_index,:)= H;
        
    end

end

function [N_values, H_values] = clustering_withdistancematrixes(distance_matrixes,ordering,z,sigtest,groupings)



q_strings = {'0','p05','p1','p2','p3','p5','1','2','5','10','20'};

N_values = cell(1, length(q_strings));
H_values = zeros(1, length(q_strings));

number_of_categories = length(groupings)/2;
n = groupings(end) - groupings(1) +1;

for q_index = 1:length(q_strings)

        distance_matrix_unordered = distance_matrixes{q_index};
        
        number_of_trials = length(distance_matrix_unordered);
        
        
        if(sigtest==1)
            ordering = randperm(number_of_trials);
        end
        
        
        distance_matrix = zeros(number_of_trials,number_of_trials);
        for i = 1:number_of_trials
            for j = i+1:number_of_trials
                distance_matrix(i,j) = distance_matrix_unordered(ordering(i),ordering(j));
                distance_matrix(j,i) = distance_matrix_unordered(ordering(i),ordering(j));

            end
        end
        
         
        
        costs = zeros(n, number_of_categories);


        for i = 1:n

            sums = zeros(1,number_of_categories);
            counts = zeros(1,number_of_categories); 

            for cat = 1:number_of_categories
                
                for j = groupings(2*cat-1):groupings(2*cat)
                    
                    if(i<j)
                          if(distance_matrix(i,j)~=0)
                        sums(cat) = sums(cat) + (distance_matrix(i,j))^(z); 
                          end
                        counts(cat) = counts(cat) + 1; 
                    else if (i>j) 
                              if(distance_matrix(j,i)~=0)
                            sums(cat) = sums(cat) + (distance_matrix(j,i))^(z);
                              end
                            counts(cat) = counts(cat)+1;
                        end
                    end  
                end
            end
            
            
            sums_nonzero_mask = (sums~=0);
            cost = zeros(1,number_of_categories);
            cost(sums_nonzero_mask) = (sums(sums_nonzero_mask)./counts(sums_nonzero_mask)).^(1/z);

            costs(i,:) = cost;
        end
        
    
     

        N = getN(costs,groupings);
        H = getH(N);

        N_values{q_index}=N;
        H_values(q_index)=H;
        
end
end


function [N] = getN(costs,groupings)

% used in clustering_withdistancematrixes to get N 
number_of_categories = length(groupings)/2;

N = zeros(number_of_categories,number_of_categories);

for cat = 1:number_of_categories
    for i = groupings(2*cat-1):groupings(2*cat)
        [x,y] = min(costs(i,:));
        count = sum(costs(i,:)==x);
        
        N(cat,:)=N(cat,:)+(costs(i,:)==x)./count;
    end
end



end




function [H] = getH(N)

[n1,n2]=size(N);

one = sum(N,1);
two = sum(N,2);

total =0;
for i=1:n1
    for j=1:n2
       
        if(N(i,j)~=0)
        total = total+ N(i,j)*(log2(N(i,j))+log2(sum(one))-log2(one(j))-log2(two(i)));
        end
    end
end

H = total/sum(one);
end


function  H_allzs_randoms = sig_test_random(n_random,distance_matrixes,groupings)

zs = [-8,-4,-2,-1,1,2,4,8];

    H_allzs_randoms = cell(n_random,1);
    
    for trial = 1:n_random

        status = trial

        H_allzs_random = zeros(8,11);
        
        for z = 1:8
            [N,H] = clustering_withdistancematrixes(distance_matrixes,0,zs(z),1,groupings);
            H_allzs_random(z,:)=H;
        end
        
        H_allzs_randoms{trial} = H_allzs_random;

    end
    


end

function [H_random_by_zs] = separate_H_randoms(H_randoms)

    
    n_random = size(H_randoms,1);
    
    for i = 1:n_random;
        claire = H_randoms{i};
        
            H_random_zm8(i,:) = claire(1,:);
            H_random_zm4(i,:) = claire(2,:);
            H_random_zm2(i,:) = claire(3,:);
            H_random_zm1(i,:) = claire(4,:);

            H_random_zp1(i,:) = claire(5,:);
            H_random_zp2(i,:) = claire(6,:);
            H_random_zp4(i,:) = claire(7,:);
            H_random_zp8(i,:) = claire(8,:);
           
    end
    
    H_random_by_zs = {H_random_zm8;
                        H_random_zm4;
                        H_random_zm2;
                        H_random_zm1;
                        H_random_zp1;
                        H_random_zp2;
                        H_random_zp4;
                        H_random_zp8;};
                    
                       
    
end



    