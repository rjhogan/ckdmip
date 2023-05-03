v = [16 32];
plot_type = '';

if 0
for iv = 1:length(v)
  model_setting{iv} = ['ecCKD-RGB-' num2str(v(iv)) ''];
end
evaluate_ckd2('ecckd-1.2',model_setting,'sw','climate','rgb-16',plot_type,'','rgb-32');
end

if 0
for iv = 1:length(v)
  model_setting{iv} = ['ecCKD-RGB-' num2str(v(iv)) '-v1'];
end
evaluate_ckd2('ecckd-1.2',model_setting,'sw','climate','rgb-16-scaled',plot_type,'','rgb-32-scaled');
else
for iv = 1:length(v)
  model_setting{iv} = ['ecCKD-RGB-' num2str(v(iv)) '-v0'];
end
evaluate_ckd2('ecckd-1.2',model_setting,'sw','climate','rgb-16-raw',plot_type,'','rgb-32-raw');
end

if 0
for iv = 1:length(v)
  model_setting{iv} = ['ecCKD-FSCK-' num2str(v(iv)) ''];
end
evaluate_ckd2('ecckd-1.2',model_setting,'lw','climate','fsck-16',plot_type,'','fsck-32');
elseif 0
for iv = 1:length(v)
  model_setting{iv} = ['ecCKD-FSCK-' num2str(v(iv)) '-v0'];
end
evaluate_ckd2('ecckd-1.2',model_setting,'lw','climate','fsck-16-raw',plot_type,'','fsck-32-raw');
end
