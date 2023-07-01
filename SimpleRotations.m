% ---- Comparing QR from householder reflections and simple rotations ---

matrixsizes = { [10, 5], [100, 50], [300,70] };
numrand = size(matrixsizes,2);
 

cyclematrices = {10,50,100};
numtest = size(cyclematrices,2);

numiterations = 100;

% What tests to run: 
% [0 0 0] are no tests
% [1 1 1 ] are all test
% first element is orthogonality test
% second element is columnspace test
% third is projection test
testarray = [1 1 1];
[orth_result, col_result, proj_result] = householdervsimpleQR(matrixsizes, numiterations, testarray, cyclematrices);

plotting(orth_result, col_result, proj_result, numrand, numtest);
%OBS,if changing the dimensions and types of testmatrices, you have to
%change the labels on the axis  with the function "xticklabels"


%The output is on the following form 

%For the orthogonaltiy test
%orth_result = [ mean of the error for householder, mean of the error for
%simple rotations, variance of the error for householder, varinace for the
%error for simple rotations] 
%Each row is for different size of random matrix
%For the predefined test matrices, only the value is given

%For the column space test
%col_result = [ mean of the error for householder, mean of the error for
%simple rotations, variance of the error for householder, varinace for the
%error for simple rotations] 
%Each row is for different size of random matrix
%For the predefined test matrices, only the value is given


%For the projection test
%proj_result = [ [mean of the error for householder], [mean of the error for
%simple rotations], [variance of the error for householder], [varinace for the
%error for simple rotations]]
%each element is a vector where the first element is for the test QQ'A and
%the second for (I-QQ')A
%Each row is for different size of random matrix
%For the predefined test matrices, only the value is given


function [orth_result, col_result, proj_result] = householdervsimpleQR(matrixsizes, numiterations, testarray, testmatrices)%add testmatrices as input

    
    numsizes = numel(matrixsizes)+numel(testmatrices);
    
    orth_result = cell(numsizes, 4);
    col_result = cell(numsizes, 4);
    proj_result = cell(numsizes, 4);
    
    
    for i=1:numsizes %do the test for all matrix sizes
        orthogonal_error_house = zeros(numiterations,1);
        orthogonal_error_simple = zeros(numiterations,1);
        
        column_error_house = zeros(numiterations,1);
        column_error_simple = zeros(numiterations,1);
        
        proj_error_house = zeros(numiterations,2);
        proj_error_simple = zeros(numiterations,2);
        
        for j=1:numiterations %number of tests per size
            if i <= numel(matrixsizes)
                 A = randn(matrixsizes{i});
            else
                A = gallery("cycol",testmatrices{i-numel(matrixsizes)});
            end
           
                [Q_house, R_house] = qr(A);
                [Q_simple, R_simple] = SimpleQR(A);
                
                %ORTHOGONALITY TEST
                if testarray(1)
                    orthogonal_error_house(j) = norm(eye(size(A,1))-Q_house*Q_house',2);
                    orthogonal_error_simple(j) = norm(eye(size(A,1))-Q_simple*Q_simple',2);
                end
                
                %COLOMN SPACE
                if testarray(2)
                    
                    colA = orth((A));

                    colhouse = orth((Q_house));
                    colsimple = orth((Q_simple));

                    column_error_house(j)  = subspace(colA,colhouse);
                    column_error_simple(j) = subspace(colA,colsimple);
                end
                
                %PROJECTION
                if testarray(3)
                    disp(i);
                    proj_error_house(j,1) = norm(A - Q_house*inv(Q_house'*Q_house)*Q_house'* A,2);
                    proj_error_simple(j,1) = norm(A - Q_simple*inv(Q_simple'*Q_simple)*Q_simple'* A,2);
                    proj_error_house(j,2) = norm(0 - (eye(size(A,1)) - Q_house'*inv(Q_house*Q_house')*Q_house)* A,2);
                    proj_error_simple(j,2) = norm(0 - (eye(size(A,1)) -Q_simple'*inv(Q_simple*Q_simple')*Q_simple)* A,2);
                end
                
            
        end
        %Create mean and variance and add it to cell array
        if testarray(1)
                orth_result{i,1} = mean(orthogonal_error_house);
                orth_result{i,2} = mean(orthogonal_error_simple);
                orth_result{i,3} = var(orthogonal_error_house);
                orth_result{i,4} = var(orthogonal_error_simple); 
               

        end
        
        if testarray(2)
                col_result{i,1} = mean(column_error_house);
                col_result{i,2} = mean(column_error_simple);
                col_result{i,3} = var(column_error_house);
                col_result{i,4} = var(column_error_simple); 

        end
        
        if testarray(3)
                proj_result{i,1} = mean(proj_error_house);
                proj_result{i,2} = mean(proj_error_simple);
                proj_result{i,3} = var(proj_error_house);
                proj_result{i,4} = var(proj_error_simple);  

        end
        

        
    end
    
end


function [Q, R] = SimpleQR(A)
  [N, M] = size(A);

    A_aug = A;

    S_cell = cell(N-1,1); 

    Q = eye(N);
    R = A;

    for i=1:M-1
        A_i = A_aug(i:end, i:end);
        x = A_i(:,1);

        length = size(x,1);
        x = x / norm(x);
        x_1 = x(1);
        x_2 = x(2:end);

        K = eye(length-1) - (x_2*x_2'/(1+x_1));
        S_x = [x_1 , x_2' ; -x_2, K];

        %create the full S
        if i == 1
            S = S_x;
        else
            S = [eye(i-1), zeros(i-1,N-(i-1)); zeros(N-(i-1),i-1), S_x];
        end

        S_cell{i} = S;

        A_aug = S *A_aug;
        Q = S*Q;
        R = S*R;

    end
  
end

function[] = plotting(orth_result, col_result, proj_result, numrand, numtest)
    close all;
    % Define the x-values
    xrand = linspace(1, numrand, numrand);
    xtest = linspace(1, numtest, numtest);

    %ORTHOGONALTIY - RANDOM
    % Define the y-values for Householder
    y1 = cell2mat( orth_result(1:numrand,1));

    % Define the variances for Householder
    var1 = cell2mat( orth_result(1:numrand,3));

    % Define the y-values for simple
    y2 = cell2mat( orth_result(1:numrand,2));

    % Define the variances for simple
    var2 = cell2mat( orth_result(1:numrand,4));

    % Create the figure and axes
    figure;
    hold on;

    errorbar( xrand,y1, var1.^(1/2), '-bo', 'LineWidth', 1.5);

    errorbar( xrand,y2, var2.^(1/2), '-ro', 'LineWidth', 1.5);

    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Orthogonality for random matrices');
    xticks(xrand)
    xticklabels({'10x5', '100x50', '300x70'});
    legend( 'Householder Reflections', 'Simple Rotations');

    grid on;

    hold off;
    f = gcf;
    exportgraphics(f,'plot1.png','Resolution',300)

    
    % Define the y-values for Householder
    y1 = cell2mat( orth_result(numrand+1:end,1));

    % Define the variances for Householder
    var1 = cell2mat( orth_result(numrand+1:end,3));

    % Define the y-values for simple
    y2 = cell2mat( orth_result(numrand+1:end,2));

    % Define the variances for simple
    var2 = cell2mat( orth_result(numrand+1:end,4));

    figure;
    hold on;


    errorbar( xtest,y1, var1.^(1/2), '-bo', 'LineWidth', 1.5);
    errorbar( xtest,y2, var2.^(1/2), '-ro', 'LineWidth', 1.5);

    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Orthogonality for Cyclic matrices');
    xticks(xtest)
    xticklabels({'10x10', '50x50', '100x100'});
    legend('Householder Reflections', 'Simple Rotations');

    grid on;

    hold off;
    f = gcf;
    exportgraphics(f,'plot2.png','Resolution',300)

    
    
    %COLOMNSPACE RANDOM
    % Define the y-values for Householder
    y1 =cell2mat(  col_result(1:numrand,1));

    % Define the variances for Householder
    var1 = cell2mat( col_result(1:numrand,3));

    % Define the y-values for Simple
    y2 = cell2mat( col_result(1:numrand,2));

    % Define the variances for Simple
    var2 = cell2mat( col_result(1:numrand,4));

    figure;
    hold on;

    errorbar( xrand,y1, var1.^(1/2), '-bo', 'LineWidth', 1.5);
    errorbar( xrand,y2, var2.^(1/2), '-ro', 'LineWidth', 1.5);

    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Columns space for random matrices');
    xticks(xrand)
    xticklabels({'10x5', '100x50', '300x70'});

    legend('Householder Reflections', 'Simple Rotations');
    grid on;

    hold off;
    f = gcf;
    exportgraphics(f,'plot3.png','Resolution',300)

    
    %Cyclic
    

    % Define the y-values for Householder
    y1 =cell2mat(  col_result(numrand+1:end,1));

    % Define the variances for Householder
    var1 = cell2mat( col_result(numrand+1:end,3));

    % Define the y-values for Simple
    y2 = cell2mat( col_result(numrand+1:end,2));

    % Define the variances for Simple
    var2 = cell2mat( col_result(numrand+1:end,4));
    
    

    figure;
    hold on;

    errorbar( xtest,y1, var1.^(1/2), '-bo', 'LineWidth', 1.5);
    errorbar( xtest,y2, var2.^(1/2), '-ro', 'LineWidth', 1.5);

    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Column space for Cyclic matrices');
    xticks(xtest)
    xticklabels({'10x10', '50x50', '100x100'});

    legend('Householder Reflections', 'Simple Rotations');

    grid on;

    hold off;
    f = gcf;
    exportgraphics(f,'plot4.png','Resolution',300)

    
    
    
    %PROJECTION MATRIX
    
    
    % Define the y-values for Householder
    y1 =cell2mat(  proj_result(1:numrand,1));
    y1_1 = y1(:,1);
    y1_2 = y1(:,2);

    % Define the variances for Householder
    var1 = cell2mat( proj_result(1:numrand,3));
    var1_1 = var1(:,1);
    var1_2 = var1(:,2);

    % Define the y-values for Simple
    y2 = cell2mat( proj_result(1:numrand,2));
    y2_1 = y2(:,1);
    y2_2 = y2(:,2);

    % Define the variances for Simple
    var2 = cell2mat( proj_result(1:numrand,4));
    var2_1 = var2(:,1);
    var2_2 = var2(:,2);

    figure;
    hold on;

    errorbar( xrand,y1_1, var1_1.^(1/2), '-bo', 'LineWidth', 1.5);
    errorbar( xrand,y2_1, var2_1.^(1/2), '-ro', 'LineWidth', 1.5);

    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Projection onto col(Q) for randn mtx');
    xticks(xrand)
    xticklabels({'10x5', '100x50', '300x70'});

    legend('Householder Reflections', 'Simple Rotations');
    grid on;

    hold off;
    f = gcf;
    exportgraphics(f,'plot5.png','Resolution',300)

    
    figure;
    hold on;

    errorbar( xrand,y1_2, var1_2.^(1/2), '-bo', 'LineWidth', 1.5);
    errorbar( xrand,y2_2, var2_2.^(1/2), '-ro', 'LineWidth', 1.5);

    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Projection onto null(Q) for randn mtx');
    xticks(xrand)
    xticklabels({'10x5', '100x50', '300x70'});

    legend('Householder Reflections', 'Simple Rotations');
    grid on;

    hold off;
    f = gcf;
    exportgraphics(f,'plot6.png','Resolution',300)

    
    %Cyclic
    
    % Define the y-values for Householder
    y1 = cell2mat( proj_result(numrand+1:end,1));
    y1_1 = y1(:,1);
    y1_2 = y1(:,2);
    
    var1 = cell2mat( proj_result(numrand+1:end,3));
    var1_1 = var1(:,1);
    var1_2 = var1(:,2);
    
    % Define the y-values for Simple
    y2 = cell2mat( proj_result(numrand+1:end,2));
    y2_1 = y2(:,1);
    y2_2 = y2(:,2);
    
    var2 = cell2mat( proj_result(numrand+1:end,4));
    var2_1 = var2(:,1);
    var2_2 = var2(:,2);
    
    figure;
    hold on;

    errorbar( xrand,y1_1, var1_1.^(1/2), '-bo', 'LineWidth', 1.5);
    errorbar( xrand,y2_1, var2_1.^(1/2), '-ro', 'LineWidth', 1.5);

    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Projection onto col(Q) for Cyclic mtx');
    xticks(xtest)
    xticklabels({'10x10', '50x50', '100x100'});

    legend('Householder Reflections', 'Simple Rotations');

    grid on;

    hold off;
    f = gcf;
    exportgraphics(f,'plot7.png','Resolution',300)

    
    figure;
    hold on;

    errorbar( xrand,y1_2, var1_2.^(1/2), '-bo', 'LineWidth', 1.5);

    errorbar( xrand,y2_2, var2_2.^(1/2), '-ro', 'LineWidth', 1.5);

    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Projection onto null(Q) for Cyclic mtx');
    xticks(xtest)
    xticklabels({'10x10', '50x50', '100x100'});

    legend('Householder Reflections', 'Simple Rotations');

    grid on;

    hold off;
    f = gcf;
    exportgraphics(f,'plot8.png','Resolution',300)

end