function norm_floor(x, min=0.01)
	x = x / sum(x)
	for (i, el) in enumerate(x)
		# if el is too small, shave off the largest element and add to el
		if el < min
			x[argmax(x)] -= min - el
			x[i] = min
		end
	end
	return x
end