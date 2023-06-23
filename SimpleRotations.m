% ---- Comparing QR from householder reflections and simple rotations ---

matrixsizes = { [10, 5], [100, 50], [300,70] };
 
%todo - add extra test matrices as a cell
testmatrices = {hilb(10), hilb(50), hilb(100)};

numiterations = 10;

% What tests to run: 
% [0 0 0] are no tests
% [1 1 1 ] are all test
% first element is orthogonality test
% second element is columnspace test
% third is projection test
testarray = [1 1 1];
[orth_result, col_result, proj_result] = householdervsimpleQR(matrixsizes, numiterations, testarray, testmatrices);
plotting(orth_result, col_result, proj_result);
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
    if  ~isequal(size(testarray), [1 3]) 
        error('Invalid input');
    end
    %add test more tests on input 
    
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
                    %Do what is already here
                 A = randn(matrixsizes{i});
            else
                if j>1
                    break;
                end
                A = testmatrices{i-numel(matrixsizes)};
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
                    colhouse = orth((Q_house*R_house));
                    colsimple = orth((Q_simple*R_simple));

                    column_error_house(j) = norm(colA - colhouse,2);
                    column_error_simple(j) = norm(colA - colsimple,2);
                end
                
                %PROJECTION
                if testarray(3)
                    proj_error_house(j,1) = norm(A - Q_house*Q_house'* A,2);
                    proj_error_simple(j,1) = norm(A - Q_simple*Q_simple'* A,2);
                    proj_error_house(j,2) = norm(0 - (eye(size(A,1)) - Q_house*Q_house')* A,2);
                    proj_error_simple(j,2) = norm(0 - (eye(size(A,1)) - Q_simple*Q_simple')* A,2);
                end
                
            
        end
        %Create mean and variance and add it to cell array
        if testarray(1)
            if i <= numel(matrixsizes)
                orth_result{i,1} = mean(orthogonal_error_house);
                orth_result{i,2} = mean(orthogonal_error_simple);
                orth_result{i,3} = var(orthogonal_error_house);
                orth_result{i,4} = var(orthogonal_error_simple); 
            else 
                orth_result{i,1} = orthogonal_error_house(1);
                orth_result{i,2} = orthogonal_error_simple(1);
            end
        end
        
        if testarray(2)
            if i <= numel(matrixsizes)
                col_result{i,1} = mean(column_error_house);
                col_result{i,2} = mean(column_error_simple);
                col_result{i,3} = var(column_error_house);
                col_result{i,4} = var(column_error_simple); 
            else
                col_result{i,1} = column_error_house(1);
                col_result{i,2} = column_error_simple(1);
            end
        end
        
        if testarray(3)
            if i <= numel(matrixsizes)
                proj_result{i,1} = mean(proj_error_house);
                proj_result{i,2} = mean(proj_error_simple);
                proj_result{i,3} = var(proj_error_house);
                proj_result{i,4} = var(proj_error_simple);  
            else
                proj_result{i,1} = proj_error_house(1,:);
                proj_result{i,2} = proj_error_simple(1,:);
            end
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

function[] = plotting(orth_result, col_result, proj_result)
    close all;
    % Define the x-values
    x = [1,2,3];

    %ORTHOGONALTIY - RANDOM
    % Define the y-values for dataset one
    y1 = cell2mat( orth_result(1:3,1));

    % Define the variances for dataset one
    var1 = cell2mat( orth_result(1:3,3));

    % Define the y-values for dataset two
    y2 = cell2mat( orth_result(1:3,2));

    % Define the variances for dataset two
    var2 = cell2mat( orth_result(1:3,4));

    % Create the figure and axes
    figure;
    hold on;

    % Plot the first dataset with error bars
    errorbar( x,y1, var1.^(1/2), '-bo', 'LineWidth', 1.5);
    %plot(x, y1, 'b-', 'LineWidth', 1.5);

    % Plot the second dataset with error bars
    errorbar( x,y2, var2.^(1/2), '-ro', 'LineWidth', 1.5);
   % plot(x, y2, 'r-', 'LineWidth', 1.5);

    % Set the axis labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Orthogonality for random matrices');
    xticks([1,2,3])
    xticklabels({'10x5', '100x50', '300x70'});

    % Set the legend
    legend( 'Householder Reflections', 'Simple Rotations');

    % Turn the grid on
    grid on;

    % Hold off from further plotting
    hold off;
    f = gcf;
    exportgraphics(f,'plot1.png','Resolution',300)

    
    %HILBERT
    
    % Define the y-values for dataset one
    y1 = cell2mat( orth_result(4:6,1));

    % Define the y-values for dataset two
    y2 = cell2mat( orth_result(4:6,2));

    % Create the figure and axes
    figure;
    hold on;

    % Plot the first dataset with error bars
    plot( x,y1,  '-bo', 'LineWidth', 1.5);

    % Plot the second dataset with error bars
    plot( x,y2, '-ro', 'LineWidth', 1.5);

    % Set the axis labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Orthogonality for Hilbert matrices');
    xticks([1,2,3])
    xticklabels({'10x10', '50x50', '100x100'});

    % Set the legend
    legend('Householder Reflections', 'Simple Rotations');

    % Turn the grid on
    grid on;

    % Hold off from further plotting
    hold off;
    f = gcf;
    exportgraphics(f,'plot2.png','Resolution',300)

    
    
    %COLOMNSPACE RANDOM
    % Define the y-values for dataset one
    y1 =cell2mat(  col_result(1:3,1));

    % Define the variances for dataset one
    var1 = cell2mat( col_result(1:3,3));

    % Define the y-values for dataset two
    y2 = cell2mat( col_result(1:3,2));

    % Define the variances for dataset two
    var2 = cell2mat( col_result(1:3,4));

    % Create the figure and axes
    figure;
    hold on;

    % Plot the first dataset with error bars
    errorbar( x,y1, var1.^(1/2), '-bo', 'LineWidth', 1.5);

    % Plot the second dataset with error bars
    errorbar( x,y2, var2.^(1/2), '-ro', 'LineWidth', 1.5);

    % Set the axis labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Columns space for random matrices');
    xticks([1,2,3])
    xticklabels({'10x5', '100x50', '300x70'});

    % Set the legend
    legend('Householder Reflections', 'Simple Rotations');
    % Turn the grid on
    grid on;

    % Hold off from further plotting
    hold off;
    f = gcf;
    exportgraphics(f,'plot3.png','Resolution',300)

    
    %HILBERT
    
    % Define the y-values for dataset one
    y1 = cell2mat( col_result(4:6,1));

    % Define the y-values for dataset two
    y2 = cell2mat( col_result(4:6,2));

    % Create the figure and axes
    figure;
    hold on;

    % Plot the first dataset with error bars
    plot( x,y1,  '-bo', 'LineWidth', 1.5);

    % Plot the second dataset with error bars
    plot( x,y2, '-ro', 'LineWidth', 1.5);

    % Set the axis labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Column space for Hilbert matrices');
    xticks([1,2,3])
    xticklabels({'10x10', '50x50', '100x100'});

    % Set the legend
    legend('Householder Reflections', 'Simple Rotations');

    % Turn the grid on
    grid on;

    % Hold off from further plotting
    hold off;
    f = gcf;
    exportgraphics(f,'plot4.png','Resolution',300)

    
    
    
    %PROJECTION MATRIX
    
    
    % Define the y-values for dataset one
    y1 =cell2mat(  proj_result(1:3,1));
    y1_1 = y1(:,1);
    y1_2 = y1(:,2);

    % Define the variances for dataset one
    var1 = cell2mat( proj_result(1:3,3));
    var1_1 = var1(:,1);
    var1_2 = var1(:,2);

    % Define the y-values for dataset two
    y2 = cell2mat( proj_result(1:3,2));
    y2_1 = y2(:,1);
    y2_2 = y2(:,2);

    % Define the variances for dataset two
    var2 = cell2mat( proj_result(1:3,4));
    var2_1 = var2(:,1);
    var2_2 = var2(:,2);

    % Create the figure and axes
    figure;
    hold on;

    % Plot the first dataset with error bars
    errorbar( x,y1_1, var1_1.^(1/2), '-bo', 'LineWidth', 1.5);

    % Plot the second dataset with error bars
    errorbar( x,y2_1, var2_1.^(1/2), '-ro', 'LineWidth', 1.5);

    % Set the axis labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Projection with QQ^T for random mtx');
    xticks([1,2,3])
    xticklabels({'10x5', '100x50', '300x70'});

    % Set the legend
    legend('Householder Reflections', 'Simple Rotations');
    % Turn the grid on
    grid on;

    % Hold off from further plotting
    hold off;
    f = gcf;
    exportgraphics(f,'plot5.png','Resolution',300)

    
    % Create the figure and axes
    figure;
    hold on;

    % Plot the first dataset with error bars
    errorbar( x,y1_2, var1_2.^(1/2), '-bo', 'LineWidth', 1.5);

    % Plot the second dataset with error bars
    errorbar( x,y2_2, var2_2.^(1/2), '-ro', 'LineWidth', 1.5);

    % Set the axis labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Projection with I-QQ^T for random mtx');
    xticks([1,2,3])
    xticklabels({'10x5', '100x50', '300x70'});

    % Set the legend
    legend('Householder Reflections', 'Simple Rotations');
    % Turn the grid on
    grid on;

    % Hold off from further plotting
    hold off;
    f = gcf;
    exportgraphics(f,'plot6.png','Resolution',300)

    
    %HILBERT
    
    % Define the y-values for dataset one
    y1 = cell2mat( proj_result(4:6,1));
    y1_1 = y1(:,1);
    y1_2 = y1(:,2);
    
    % Define the y-values for dataset two
    y2 = cell2mat( proj_result(4:6,2));
    y2_1 = y2(:,1);
    y2_2 = y2(:,2);
    
    % Create the figure and axes
    figure;
    hold on;

    % Plot the first dataset with error bars
    plot( x,y1_1,  '-bo', 'LineWidth', 1.5);

    % Plot the second dataset with error bars
    plot( x,y2_1, '-ro', 'LineWidth', 1.5);

    % Set the axis labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Projection with QQ^T for Hilbert mtx');
    xticks([1,2,3])
    xticklabels({'10x10', '50x50', '100x100'});

    % Set the legend
    legend('Householder Reflections', 'Simple Rotations');

    % Turn the grid on
    grid on;

    % Hold off from further plotting
    hold off;
    f = gcf;
    exportgraphics(f,'plot7.png','Resolution',300)

    
    % Create the figure and axes
    figure;
    hold on;

    % Plot the first dataset with error bars
    plot( x,y1_2,  '-bo', 'LineWidth', 1.5);

    % Plot the second dataset with error bars
    plot( x,y2_2, '-ro', 'LineWidth', 1.5);

    % Set the axis labels and title
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Error in Projection of with I-QQ^T for Hilbert mtx');
    xticks([1,2,3])
    xticklabels({'10x10', '50x50', '100x100'});

    % Set the legend
    legend('Householder Reflections', 'Simple Rotations');

    % Turn the grid on
    grid on;

    % Hold off from further plotting
    hold off;
    f = gcf;
    exportgraphics(f,'plot8.png','Resolution',300)

    
    

end
%https://laurenthoeltgen.name/post/qr-benchmark/