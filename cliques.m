function [L,C,S,R,par] = cliques(A, varargin)
% CLIQUES   Computes cliques and clique tree of filled sparse matrix.
%    [L,C,S,R,par] = CLIQUES(A) returns a sparse lower triangular
%    matrix L with the filled pattern, a cell array C with the cliques, a
%    cell array S with the separator sets, a cell array R with the
%    residual sets, and a parent array par. The input matrix A must
%    be a sparse symmetric matrix.
%
%    [L,C,S,R,par] = CLIQUES(A,p) is the same as [L,C,S,R,par] =
%    CLIQUES(A(p,p)). 
%   
%    The parent array par has the following meaning: if par(i) = 0, 
%    then clique i is a root clique in the clique forest, and
%    otherwise par(i) is the index of the parent of clique i.

n = size(A,1);
assert(n == size(A,2),'A must be square')
if nargin  == 1
    p = 1:n;
elseif nargin == 2
    p = varargin{1};
else
    error('Wrong number of arguments.')
end


% compute symbolic factorization and run Pothen-Sun alg.
[colcount, h, parent, post, L] = symbfact(A(p,p),'sym','lower');
[snodes, flag] = pothen_sun(parent, post, colcount);

% extract representative vertices
rv = find(flag < 0);

% extract supernodes (residual sets)
R = cell(n,1);
Rn = zeros(n,1);  % array for counting number of nodes added to supernodes 
for ii = 1:n
    f = flag(ii);
    if f < 0
        % initialize array if vertex ii is a representative vertex
        R{ii} = zeros(-f,1); 
        R{ii}(1) = ii;
        Rn(ii) = Rn(ii)+1;
    end
end
for ii = 1:n
    f = flag(ii);
    if f > 0
        % vertex ii belongs to supernode with representative vertex f
        Rn(f) = Rn(f) + 1;  
        R{f}(Rn(f)) = ii;
    end
end

% extract cliques
C = cell(n,1);
for ii = 1:n
    f = flag(ii);
    if f < 0
        C{ii} = find(L(:,ii) ~= 0);
    end
end

% extract separator sets
S = cell(n,1);
for ii = 1:n
    f = flag(ii);
    if f < 0
        S{ii} = setdiff(C{ii}, R{ii});
    end
end

% find supernodal parent structure
snpar = zeros(n,1);
for ii = 1:snodes
    v = rv(ii);
    pi = parent(v); % parent of node v in etree
    while pi > 0                 
        if flag(pi) ~= v    % parent belongs to different supernode
            if flag(pi) > 0
                snpar(v) = flag(pi);
            else
                snpar(v) = pi;
            end
            break                
        else % parent is not a representative vertex
            pi = parent(pi);
        end
    end
end

% compress cell arrays 
S = S(rv); % cell array with separator sets
R = R(rv); % cell array with residual sets
C = C(rv); % cell array with cliques

% compress supernodal parent array
reverse_rv = zeros(n,1);
reverse_rv(rv) = 1:snodes;
par = zeros(snodes,1);
for ii = 1:snodes
    if snpar(rv(ii)) > 0
        par(ii) = reverse_rv(snpar(rv(ii)));
    else
        par(ii) = 0;
    end
end


function [ns, flag] = pothen_sun(par, post, colcount)
% POTHEN_SUN   Finds the representative vertices and the number of supernodes.
%    [n,flag] = POTHEN_SUN(par,post,colcount) computes the
%    number of supernodes ns and an array flag  of length n; if
%    flag(i) < 0, then -flag(i) is the degree of the supernode with
%    representative vertex i, and if flag(i) >= 0, then flag(i) is
%    the representative vertex to which node i belongs.
%            
 
n = length(par);
flag = -ones(n,1);
ns = n;
for j = post
    mdeg = colcount(j) - 1;
    if par(j) ~= 0 && mdeg == colcount(par(j)) && flag(par(j)) == -1
        % par(j) not assigned to supernode
        ns = ns - 1;
        if flag(j) < 0    % j is a repr. vertex
            flag(par(j)) = j;
            flag(j) = flag(j) -  1;
        else              % j is not a repr. vertex
            flag(par(j)) = flag(j);
            flag(flag(j)) = flag(flag(j)) - 1;    
        end    
    end
end
