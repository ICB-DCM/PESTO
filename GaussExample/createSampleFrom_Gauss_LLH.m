define_Gauss_LLH();

sample=[0,0];
for i = 1:100000
  u = rand(); 
  x=-2 +rand()*14;
  y=-2 +rand()*14; 
  if (exp(logP([x,y]'))>u)
    sample(end+1,:)=[x,y];
  end
end
sample = sample(2:end,:);

plot(sample(:,1),sample(:,2),'.');
hold all;
plot_Gauss_LLH();