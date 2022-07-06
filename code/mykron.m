%% mykron.m
% Evan Fortunato
% April 2000
%
% Description below by original author

function PO=mykron(varargin)
%function PO=mykron(varargin)
% This little function allows you to pass all of
% the product operators in at once and it goes through
% and makes the right answer.  It is equivalent to
% kron(kron(kron(x,x),x),x)....
% written April 29th 2000, -Evan Fortunato

if length(varargin) < 1
	error('Please enter atleast 1 Product Operator!');
elseif length(varargin) == 1
	PO = varargin{1};
else
	PO = varargin{1};
	for j=2:length(varargin)
		PO=kron(PO, varargin{j});
    end
end


