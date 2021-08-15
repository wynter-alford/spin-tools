function mq=mqc(N);

mq0=[0 -1;1 0];
mq=mq0;
for k=2:N
  mq=[mq mq+mq0(1,2)*ones(2^(k-1),2^(k-1));mq+mq0(2,1)*ones(2^(k-1),2^(k-1)) mq];
end
