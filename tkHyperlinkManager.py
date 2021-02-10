from tkinter import *

class HyperlinkManager:

    def __init__(self, text):

        self.text = text
        self.text.tag_config("hyper", foreground="blue", underline=1)
        self.text.tag_bind("hyper", "<Button-1>", self._click)

        self.reset()

    def reset(self):
        self.links = {}
        self.url = {}

    def add(self, action, url):
        tag = "hyper-%d" % len(self.links)
        self.links[tag] = action
        self.url[tag] = url
        return "hyper", tag

    def _click(self, event):
        for tag in self.text.tag_names(CURRENT):
            if tag[:6] == "hyper-":
                self.links[tag](self.url[tag])
                return
