function C = mmat(A,B,dim)
% Simple matrix multiplication of multidimensional arrays.
%
% Input:
% A, B      Multidimensional input arrays.
% dim       Contains two number, that selects two dimensions.
%
% The multiplication is a standard matrix multiplication. The two matrices
% are selected by dim:
%   AB = A(... dim(1) ... dim(2) ...) * B(... dim(1) ... dim(2) ...)
% The necessary condition that the multiplication can be performed:
%   size(A,dim(2)) = size(B,dim(1))
%
% Singleton dimensions in both A and B matrices are supported.
%
% The default value for dim is dim = [1 2].
%
% Examples:
% For 2D matrices mmat is identical to the Matlab built-in multiplication:
% A = [1 2 3 4];
% B = [1;2;3;4];
% C = mmat(A,B)
%
% C will be 30.
%
% For multidimensional arrays:
% A = repmat([1 2 3 4],[1 1 5]);
% B = [1 2 3 4]';
% C = mmat(A,B)
% C will be an array with dimensions of 1x1x5 and every element is 30.
% Copyright (c) 2013, Sándor Tóth
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if nargin == 0
    help mmat;
    return;
end

if (nargin < 3)
    dim = [1 2];
end

if numel(dim)~=2
    error('sw:sw_matmult:WrongInput','dim has to be a two element array!');
end

if size(A,dim(2)) ~= size(B,dim(1))
    error('sw:sw_matmult:WrongInput','Wrong input matrix sizes!');
end

nDA = ndims(A);
nDB = ndims(B);
nD = max(nDA,nDB);

nA = [size(A),ones(1,nD-nDA)]; nA = nA(dim); 
nB = [size(B),ones(1,nD-nDB)]; nB = nB(dim);

% form A matrix
% (nA1) x (nA2) x nB2
A = repmat(A,[ones(1,nD) nB(2)]);
% form B matrix
% nA1 x (nB1) x (nB2)
idx = 1:nD+1; idx([dim end]) = idx([end dim]);
repv = ones(1,nD+1); repv(dim(1)) = nA(1);

B = repmat(permute(B,idx),repv);

% multiply with expanding along singleton dimensions
C = sum(bsxfun(@times,A,B),dim(2));

idx2 = 1:nD+1; idx2([dim end]) = idx2([dim(1) end dim(2)]);

% permute back the final result to the right size
C = permute(C,idx2);

end