for i=1:10
    A = gallery("cycol",100);
    [q, r] = SimpleQR(A);
    hei = orth(q);

end


function [Q, R] = SimpleQR(A)
  [N, M] = size(A);

    A_aug = A;

    S_cell = cell(N-1,1); 

    Q = eye(N);
    R = A;

    for i=1:M-1
        A_i = A_aug(i:end, i:end);
        if A_i(:,1) == zeros(size(A_i(:,1)))
            break;
        end
        x = A_i(:,1);
        if isnan( x )
            disp("halla");
        end
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