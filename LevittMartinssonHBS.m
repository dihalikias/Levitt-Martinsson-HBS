function [A, U, V, D, D0] = LevittMartinssonHBS(Ax, N, r)
% Learn a symmetric rank-r HBS matrix with matrix-vector products
%
% A = LevittMartinssonHBS(Ax, N, r) learns a symmetric rank-r NxN HBS matrix via matrix-vector
% products. Here, AX is a function handle that performs a matrix-vector
% product involving A. 
%
% The algorithm is based on the HBS algorithm that can be found in this
% paper: https://arxiv.org/abs/2205.02990
%

p = 10; % oversampling parameter

% Martinsson says to take s = 3 * r', where r' = r + p, to ensure that taking the
% pseudoinverses of matrices is stable.
s = 3 * (r + p);

% How many levels of partitioning are there?
n = log2( N );

% For now, set the level L to be as large as possible, i.e. the finest
% possible level of partitioning.
L = n - floor(log2(r)) - 1;
m = 2^(n - L);

% Make a zero  matrix to store the current guess:
A = zeros(N, N);

% Get s random samples of A and A'
Omega = randn(N, s);
Y = Ax(Omega);

% Make cell arrays for each level's diagonal block matrices U, V, and D,
% where the matrix U at a given level is U{level}
U = cell(L, 1);
V = cell(L, 1);
D = cell(L, 1);

% Base case in telescoping factorization
D0 = zeros(2 * r, 2 * r);

for level = L:-1:0
    LastLevel = level + 1;
            
    % Index the submatrices of U(level), V(level), and D(level) corresponding to this level
    if level == L
        I1 = reshape(1:N, m, N/m);
        UL = zeros(N, 2^L * r);
        VL = zeros(N, 2^L * r);
        DL = zeros(N, N);   
    else
        Y = U{LastLevel}' * (Y - (D{LastLevel} * Omega));
            
        Omega = V{level+1}' * Omega;
        
        I1 = reshape(1: r * 2^(level + 1), 2*r, 2^(level));
        
        % Create new block basis matrices for this level
        UL = zeros(2 * r * 2^(level), r * 2^(level));
        VL = zeros(2 * r * 2^(level), r * 2^(level));
        DL = zeros(2 * r * 2^(level), 2 * r * 2^(level));

    end
    
    % Index the columns for the previous level 
    I2 = reshape(1:r * 2^level, r, 2^(level));

    % Go through all the subblocks at a given level to form U(level) and V(level):
    for j = 1 : 2^(level)
        
        % Indices for new block matrices at this level
        subI1 = I1(:, j);
        subI2 = I2(:, j);
        
        if level == L
            OmegaT = Omega(subI1, :);

            YT = Y(subI1, :);

        else
                
            YT = Y(subI1, :);
            
            OmegaT = Omega(subI1, :);
        end

        if level > 0
            % Find nullspaces of Omega_tau and Psi_tau
            PT = null(OmegaT);
            PT = PT(:, 1:r);

            % Form an orthonormal basis for the blocks at this level by multiplying by
            % the nullspace to eliminate contributions from the diagonal
            % block
            [UT, ~, ~] = qr(YT * PT, 0);
            UT = UT(:, 1:r);
            VT = UT;
            
            % Size of identity matrix
            isize = 2 * r;
            
            if level == L
                isize = m;
            end
            
            I = eye(isize, isize);
            DT = ((I - (UT * UT')) * YT * pinv(OmegaT)) + (UT*UT'*((I - (VT * VT')) * YT * pinv(OmegaT))');
            UL(subI1, subI2) = UT;
            VL(subI1, subI2) = UT;
            
        else
            % At the root node, just return A = D(0)
            DT = YT * pinv(OmegaT);
        end

        DL(subI1, subI1) = DT;

    end

    if level > 0
         U{level} = UL;
         V{level} = VL;
         D{level} = DL;
       
    else
        D0 = DL;
    end

end

% Telescoping factorization 

A = Telescope(A, L, U, V, D, D0);

end

function A = Telescope(A, level, U, V, D, D0)

    if level > 0
        Ul = U{level};
        Vl = V{level};
        Dl = D{level};
        
        AA = Telescope(A, level - 1, U, V, D, D0);
        A = (Ul * AA * Vl') + Dl;
    else
        % Base case: tilde A(0) =  = D(0)
        A = D0; 
    end
end
