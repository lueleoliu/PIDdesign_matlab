function result = GAMultiObj(sysIn,is_pitch,is_first,first_td)

global sys base_set

%% Load Model %%
sys = sysIn;
cal_set = textread('config.txt');
gen = cal_set(1,1);
popsize = cal_set(1,2);

base_set = textread('baseline.txt');

range_set = textread('range.txt');

options = gaoptimset('Generations',gen,'PopulationSize',popsize, 'Display', 'diagnose');
% options = gaoptimset('Generations',gen,'PopulationSize',popsize);

if is_pitch
    lb = [range_set(1,1);range_set(1,2);range_set(1,3);range_set(1,4)];
    ub = [range_set(2,1);range_set(2,2);range_set(2,3);range_set(2,4)]; 
    if is_first    
        [x,f,~] = gamultiobj(@obj_p,4,[],[],[],[],lb,ub,@nonl_p,options);
    else
        Aeq = [0 0 0 1];
        beq = first_td;
        [x,f,~] = gamultiobj(@obj_p,4,[],[],Aeq,beq,lb,ub,@nonl_p,options);
    end
else
    lb = [range_set(3,1);range_set(3,2)];
    ub = [range_set(4,1);range_set(4,2)];
    [x,f,~] = gamultiobj(@obj_q,2,[],[],[],[],lb,ub,@nonl_q,options);
end

if ~is_pitch
    f_b = [100,100];
    i_b = 0;
    for i = 1:size(f,1)
        if (f(i,1) + 10*f(i,2)) <= (f_b(1) + 10*f_b(2))
            f_b(1) = f(i,1);
            f_b(2) = f(i,2);
            i_b = i;
        end
    end
    x_b = x(i_b,:);

    result = [x_b,f_b];
else
    score = zeros(1,size(f,1));
    for i = 1:size(f,1)

        if f(i,1)<15
            score(i) = score(i) + 5;
        elseif f(i,1)<20
            score(i) = score(i) + 4;
        elseif f(i,1)<25
            score(i) = score(i) + 3;
        elseif f(i,1)<30
            score(i) = score(i) + 2;
        elseif f(i,1)<35
            score(i) = score(i) + 1;
        else
            score(i) = score(i) + 0;
        end
        if f(i,2)<5
            score(i) = score(i) + 5;
        elseif f(i,2)<6
            score(i) = score(i) + 4;
        elseif f(i,2)<7
            score(i) = score(i) + 3;
        elseif f(i,2)<8
            score(i) = score(i) + 2;
        elseif f(i,2)<9
            score(i) = score(i) + 1;
        else
            score(i) = score(i) + 0;
        end
        if f(i,3)<20
            score(i) = score(i) + 5;
        elseif f(i,3)<25
            score(i) = score(i) + 4;
        elseif f(i,3)<30
            score(i) = score(i) + 3;
        elseif f(i,3)<35
            score(i) = score(i) + 2;
        elseif f(i,3)<40
            score(i) = score(i) + 1;
        else
            score(i) = score(i) + 0;
        end
        if f(i,4) == 0
            score(i) = score(i) + 0;
        elseif f(i,4)<3
            score(i) = score(i) + 5;
        elseif f(i,4)<4
            score(i) = score(i) + 4;
        elseif f(i,4)<5
            score(i) = score(i) + 3;
        elseif f(i,4)<6
            score(i) = score(i) + 2;
        elseif f(i,4)<7
            score(i) = score(i) + 1;
        else 
            score(i) = score(i) + 0;
        end
    end

    result_all = [x,f,score'];

    if exist('Result_all.txt','file')
        delete('Result_all.txt')
    end

    fid = fopen('Result_all.txt','wt');
    for i = 1:size(result_all,1)
        for j = 1:size(result_all,2)
            fprintf(fid, '%.3f ', result_all(i,j));
        end
        if i < size(result_all,1)
            fprintf(fid, '\n');
        end
    end
    fclose(fid);

    high_score = find(score == max(score));
    for i = 1:length(high_score)
        f_b = [100,100];
        j = high_score(i);
        if (f(j,1) + 10*f(j,2)) <= (f_b(1) + 10*f_b(2))
            f_b(1) = f(j,1);
            f_b(2) = f(j,2);
            j_b = j;
        end
    end

    x_b = x(j_b,:);

    result = [x_b,f_b];
end
    




