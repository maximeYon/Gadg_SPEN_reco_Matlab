function Out=OperatorTest(A,ImSize,DataSize,nChannels,Mode)
x = randn(ImSize) + 1j*randn(ImSize);
y = randn([DataSize(1:2) nChannels]) + 1j*randn([DataSize(1:2) nChannels]);
Ax = A*x;
Aty = A'*y;
Out=abs(x(:)'*Aty(:) - conj(y(:)'*Ax(:)));
if(nargin<5)
    disp(['Operator conj test: ' num2str(Out)]);
else
    Err=Out;
    Out=Out<1e-10;
    if(~Out)
%         error('Failed OperatorTest');
        disp(['Operator conj test: ' num2str(Err)]);
    end
end