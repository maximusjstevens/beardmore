function result=convstack(data,R);

if (R==1);
  result=data;
else
  result=stack(data,R)./stack(ones(size(data)),R);
end

function result=stack(data, R);

centers=1:R:size(data,2);
  result=zeros(size(data,1),length(centers));
  trap=gautrap(1:(1.4*R),.6*R,.15*R)'; trap=trap./sum(trap);
  window=ones(size(data(:,1)))*trap; 
  windowlength=size(window,2); 
  for n=1:length(centers)-1
    traces=round(max(1,centers(n)-windowlength/2):min(centers(n)+windowlength/2-1));
    result(:,n)=sum(window(:,1:length(traces)).*data(:,traces),2) ;
  end 
  result(:,end)=data(:,end);

function window=gautrap(x,width,sigma);

%window=gautrap[x,width,sigma];
range=max(x)-min(x);
spacing=mean(diff(x));
window=zeros(size(x));
% g=gaussian([min(x):spacing:(max(x)-width)]',(min(x)+max(x)-width)/2,sigma);
g=normpdf([min(x):spacing:(max(x)-width)]',(min(x)+max(x)-width)/2,sigma);
window=[g(1:floor(length(g)/2));ones(floor(width/spacing),1)*max(g);g(floor(length(g)/2)+1:end)];

