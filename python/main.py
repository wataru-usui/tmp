import numpy as np
import plotly.graph_objs as go

class Distribution():
    def __init__(self, a, b):
        self.a = a
        self.b = b
        
    def pdf(self, x):
        return 0.0
    
    def cdf(self, x):
        return 0.0
    
    def quantile(self, x):
        return 0.0

class Uniform(Distribution):
    def __init__(self, a, b):
        super().__init__(a, b)
    
    def pdf(self, x):
        return 1.0 / (self.b - self.a)
    
    def cdf(self, x):
        return (x - self.a) / (self.b - self.a)
    
    def quantile(self, x):
        return (self.b - self.a) * x + self.a

def integrand(x):
    return 3.0 * x * x

def antiderivative(x):
    return x * x * x

numSimulations = 50
numSamples = 250
a = -4.0
b = 5.0
fig = go.Figure(layout=go.Layout(title="Visualizing The Convergence of A Monte Carlo Integration", xaxis=dict(title="Number of Samples"), yaxis=dict(title="Estimate")))
dist = Uniform(a, b)
for i in range(numSimulations):
    xs = list(range(1, numSamples + 1))
    ys = list(map(lambda s, n: s / n, np.cumsum(list(map(lambda x: integrand(dist.quantile(x)) / dist.pdf(dist.quantile(x)), np.random.random(numSamples)))), xs))
    fig.add_trace(go.Scatter(x=xs, y=ys, showlegend=False))
fig.add_hline(y=antiderivative(b) - antiderivative(a))
fig.write_html("./index.html")