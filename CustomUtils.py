from IPython.display import display, clear_output
import ipywidgets as widgets
class ProgressBar():
    def __init__(self, endVal = 100, name = "Progress"):
        self.w = widgets.IntProgress(
            value=0,
            min=0,
            max=endVal,
            step=1,
            description= name + ":",
            bar_style='success', # 'success', 'info', 'warning', 'danger' or ''
            orientation='horizontal'
        )
        display(self.w)
        self.p = self.progress_bar()
        
    def progress_bar(self):
        while self.w.value < self.w.max:
            self.w.value += 1
            yield
    
    def move(self):
        try: 
            self.p.__next__()
            return True
        except StopIteration:
            return False
        """if self.w.value < self.w.max:
            self.p.__next__()
            return True
        return False"""